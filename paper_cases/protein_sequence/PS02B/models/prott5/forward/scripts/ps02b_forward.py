#!/usr/bin/env python3
"""Development-only PS02B encoding, scoring, and topology-macro aggregation."""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math, re, sys, time
from pathlib import Path
import numpy as np
import torch
from safetensors.torch import save_file
from safetensors import safe_open

BASE=Path.cwd()
ROOT=None
SRC=None
FASTA=None
PANELS={}
MODELS={'esm2':None,'prott5':None,'esmc':None}

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def gzrows(p):
 with gzip.open(p,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def read_objects():
 with (SRC/'OBJECT_CANDIDATE.tsv').open(newline='') as f:
  return {r['domain_id']:r for r in csv.DictReader(f,delimiter='\t')}
def fasta(wanted):
 out={}; key=None; chunks=[]
 with FASTA.open() as f:
  for line in f:
   if line.startswith('>'):
    if key in wanted: out[key]=''.join(chunks)
    key=line[1:].split()[0].replace('cath44_s20:','').split('|')[-1].split('/')[0]; chunks=[]
   elif key in wanted: chunks.append(line.strip())
  if key in wanted:out[key]=''.join(chunks)
 return out
def prepare():
 rel={k:gzrows(p/'SELECTED_CORE_RELATIONS.tsv.gz') for k,p in PANELS.items()}
 ids={x for rows in rel.values() for r in rows for x in (r['anchor'],r['positive'],r['negative'])}
 obj=read_objects(); seq=fasta(ids); errs=[]
 for i in sorted(ids):
  if i not in obj: errs.append([i,'OBJECT_MISSING'])
  elif i not in seq:errs.append([i,'SEQUENCE_MISSING'])
  elif sha_text(seq[i])!=obj[i]['sequence_sha256']:errs.append([i,'SEQUENCE_HASH_MISMATCH'])
 ROOT.mkdir(parents=True,exist_ok=True)
 with (ROOT/'DEVELOPMENT_ENCODING_INPUT_MANIFEST.tsv').open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['domain_id','sequence_sha256','sequence_length','topology_id','hsf_id','parent_proxy'])
  for i in sorted(ids):
   if i in obj:w.writerow([i,obj[i]['sequence_sha256'],len(seq.get(i,'')),obj[i]['topology_id'],obj[i]['hsf_id'],obj[i]['parent_proxy_pdb_chain']])
 report={'kind':'PS02B_DEVELOPMENT_ENCODING_INPUT','formal_result':False,'panels':{k:{'relation_count':len(v),'relation_sha256':sha(PANELS[k]/'SELECTED_CORE_RELATIONS.tsv.gz')}for k,v in rel.items()},'expected_endpoint_count':len(ids),'found_sequence_count':len(seq),'failure_count':len(errs),'failures':errs[:20],'source_object_ledger_sha256':sha(SRC/'OBJECT_CANDIDATE.tsv'),'source_fasta_sha256':sha(FASTA),'input_manifest_sha256':sha(ROOT/'DEVELOPMENT_ENCODING_INPUT_MANIFEST.tsv')}
 (ROOT/'PREPARE_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));
 if errs:raise SystemExit('PREPARE_BLOCKED')
 with (ROOT/'endpoint_sequences.faa').open('w') as f:
  for i in sorted(ids):f.write(f'>{i}\n{seq[i]}\n')
def sha_text(s):return hashlib.sha256(s.encode('ascii')).hexdigest()
def sequences():
 out={};key=None
 with (ROOT/'endpoint_sequences.faa').open() as f:
  for x in f:
   if x.startswith('>'):key=x[1:].strip();out[key]=''
   else:out[key]+=x.strip()
 return out
def batches(seqs,n,t):
 b=[];total=0
 for i in sorted(seqs,key=lambda x:(len(seqs[x]),x)):
  c=len(seqs[i])+2
  if b and (len(b)>=n or total+c>t):yield b;b=[];total=0
  b.append(i);total+=c
 if b:yield b
def encode(name):
 seqs=sequences();out=ROOT/name;out.mkdir(exist_ok=True);dest=out/'embeddings.safetensors'
 if dest.exists():raise SystemExit('EXISTING_OUTPUT_PRESERVED')
 torch.cuda.set_device(0);torch.cuda.reset_peak_memory_stats();torch.manual_seed(20260930)
 if name=='esm2':
  sys.path.insert(0,str(BASE/'m3_smoke/python_deps'));import esm
  m,a=esm.pretrained.load_model_and_alphabet_local(str(MODELS[name]));m=m.eval().half().cuda();conv=a.get_batch_converter()
  def fw(ids):
   _,_,tok=conv([(i,seqs[i]) for i in ids])
   with torch.inference_mode():h=m(tok.cuda(),repr_layers=[33],return_contacts=False)['representations'][33]
   return [h[j,1:len(seqs[i])+1].float().mean(0).cpu().contiguous() for j,i in enumerate(ids)]
  bi=batches(seqs,16,4096);meta={'model':'esm2_t33_650M_UR50D','layer':'33','pooling':'residue_mean'}
 elif name=='prott5':
  sys.path.insert(0,str(BASE/'m3_encoder_audit/python_deps_prott5_444'));from transformers import T5EncoderModel,T5Tokenizer
  tok=T5Tokenizer.from_pretrained(str(MODELS[name]),do_lower_case=False,legacy=True,local_files_only=True);m=T5EncoderModel.from_pretrained(str(MODELS[name]),torch_dtype=torch.float16,local_files_only=True).cuda().eval()
  def fw(ids):
   enc=tok([' '.join(re.sub(r'[UZOB]','X',seqs[i])) for i in ids],add_special_tokens=True,padding=True,return_tensors='pt')
   with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):h=m(input_ids=enc['input_ids'].cuda(),attention_mask=enc['attention_mask'].cuda()).last_hidden_state
   return [h[j,:len(seqs[i])].float().mean(0).cpu().contiguous() for j,i in enumerate(ids)]
  bi=batches(seqs,8,2048);meta={'model':'Rostlab/prot_t5_xl_half_uniref50-enc','pooling':'residue_mean','UO_mapping':'X'}
 else:
  from esm.models.esmc import EsmcForMaskedLM,EsmcTokenizer
  tok=EsmcTokenizer();m=EsmcForMaskedLM.from_pretrained(str(MODELS[name]),local_files_only=True).cuda().eval()
  def fw(ids):
   enc={k:v.cuda() for k,v in tok([seqs[i] for i in ids],return_tensors='pt',padding=True).items()}
   with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):h=m(**enc).last_hidden_state
   return [h[j,1:len(seqs[i])+1].float().mean(0).cpu().contiguous() for j,i in enumerate(ids)]
  bi=batches(seqs,16,4096);meta={'model':'biohub/ESMC-600M','pooling':'residue_mean'}
 vec={};idx=[];fail=[];start=time.time()
 for nb,ids in enumerate(bi,1):
  try:r=fw(ids)
  except Exception as e:
   r=[]
   for i in ids:
    try:r.append((i,fw([i])[0]))
    except Exception as e2:fail.append((i,type(e2).__name__,str(e2)[:180]))
   for i,x in r:
    if torch.isfinite(x).all():
     k=f'embedding_{len(vec)+1:06d}';vec[k]=x;idx.append((i,k,len(seqs[i]),x.numel(),hashlib.sha256(x.numpy().tobytes()).hexdigest()))
    else:fail.append((i,'NONFINITE_EMBEDDING',''))
   continue
  for i,x in zip(ids,r):
   if torch.isfinite(x).all():
    k=f'embedding_{len(vec)+1:06d}';vec[k]=x;idx.append((i,k,len(seqs[i]),x.numel(),hashlib.sha256(x.numpy().tobytes()).hexdigest()))
   else:fail.append((i,'NONFINITE_EMBEDDING',''))
  if nb%20==0:print(json.dumps({'model':name,'batches':nb,'encoded':len(idx),'failed':len(fail)}),flush=True)
 if vec:save_file(vec,str(dest),metadata=meta)
 with (out/'embedding_index.tsv').open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['domain_id','embedding_key','sequence_length','embedding_dim','embedding_sha256']);w.writerows(idx)
 with (out/'failure_ledger.tsv').open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['domain_id','failure_type','detail']);w.writerows(fail)
 rep={'status':'PASS' if len(idx)==len(seqs) and not fail else 'PARTIAL','kind':'PS02B_DEVELOPMENT_MODEL_FORWARD','formal_result':False,'model':name,'endpoint_count':len(seqs),'encoded_count':len(idx),'failure_count':len(fail),'checkpoint_sha256':sha(MODELS[name]) if MODELS[name].is_file() else sha(MODELS[name]/('pytorch_model.bin' if name=='prott5' else 'model.safetensors')),'embedding_file_sha256':sha(dest) if dest.exists() else None,'embedding_dim':idx[0][3] if idx else None,'seconds':time.time()-start,'peak_gpu_memory_mib':torch.cuda.max_memory_allocated()/2**20}
 (out/'FORWARD_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep),flush=True)
def cosine(a,b):
 a=np.asarray(a,dtype=np.float64);b=np.asarray(b,dtype=np.float64);ma=np.max(np.abs(a));mb=np.max(np.abs(b))
 if ma==0 or mb==0:return 0.0
 a=a/ma;a=a/np.linalg.norm(a);b=b/mb;b=b/np.linalg.norm(b);return max(-1.,min(1.,math.fsum(float(x) for x in a*b)))
def macro(rows,key):
 groups={}
 for r in rows: groups.setdefault(tuple(r[k] for k in key),[]).append(r)
 result=[]
 for g,items in groups.items():
  x=dict(items[0]); x.update({k:v for k,v in zip(key,g)}); x['value']=math.fsum(y['value'] for y in items)/len(items); x['n_child']=len(items); result.append(x)
 return result
def score(name):
 d=ROOT/name;idx={r['domain_id']:r['embedding_key'] for r in csv.DictReader((d/'embedding_index.tsv').open(),delimiter='\t')}
 with safe_open(str(d/'embeddings.safetensors'),framework='np') as f:emb={i:f.get_tensor(k) for i,k in idx.items()}
 obj=read_objects();summ={}
 for label,panel in PANELS.items():
  relations=gzrows(panel/'SELECTED_CORE_RELATIONS.tsv.gz');out=d/f'PAIR_SCORES_{label}_DEVELOPMENT.tsv.gz';fails=[];leaf=[]
  with gzip.open(out,'wt',newline='') as f:
   w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['anchor','positive','negative','topology','anchor_hsf','negative_hsf','anchor_parent','positive_parent','negative_parent','positive_cosine','negative_cosine','margin','credit'])
   for r in relations:
    try:
     ps=cosine(emb[r['anchor']],emb[r['positive']]);ns=cosine(emb[r['anchor']],emb[r['negative']]);margin=ps-ns;credit=1. if margin>0 else -1. if margin<0 else 0.
     x={'anchor':r['anchor'],'positive':r['positive'],'negative':r['negative'],'topology':r['topology'],'anchor_hsf':r['anchor_hsf'],'negative_hsf':r['negative_hsf'],'anchor_parent':obj[r['anchor']]['parent_proxy_pdb_chain'],'positive_parent':obj[r['positive']]['parent_proxy_pdb_chain'],'negative_parent':r['negative_parent'],'value':credit}
     leaf.append(x);w.writerow([x[k] for k in ('anchor','positive','negative','topology','anchor_hsf','negative_hsf','anchor_parent','positive_parent','negative_parent')]+[repr(ps),repr(ns),repr(margin),repr(credit)])
    except Exception as e:fails.append([r['anchor'],r['negative'],str(e)[:150]])
  a=macro(leaf,['anchor','positive','negative_hsf','negative_parent'])
  b=macro(a,['anchor','positive','negative_hsf','negative_parent'])
  c=macro(b,['anchor','positive','negative_hsf'])
  d1=macro(c,['anchor','positive'])
  e=macro(d1,['anchor','positive_parent'])
  f1=macro(e,['anchor','positive_parent'])
  g=macro(f1,['anchor'])
  h=macro(g,['topology','anchor_hsf','anchor_parent'])
  i=macro(h,['topology','anchor_hsf','anchor_parent'])
  j=macro(i,['topology','anchor_hsf'])
  k=macro(j,['topology'])
  overall=macro(k,[])[0]['value'] if k else None
  with (d/f'TOPOLOGY_SCORES_{label}_DEVELOPMENT.tsv').open('w',newline='') as q:
   w=csv.writer(q,delimiter='\t',lineterminator='\n');w.writerow(['topology','F','anchor_hsf_count']);w.writerows([(x['topology'],repr(x['value']),x['n_child']) for x in k])
  summ[label]={'relation_count':len(relations),'pair_scores_written':len(leaf),'failure_count':len(fails),'pair_score_sha256':sha(out),'F_PS02B':overall,'topology_count':len(k),'anchor_hsf_count':len(j),'anchor_source_unit_count':len(i),'anchor_count':len(g),'negative_hsf_cell_count':len(c)}
  (d/f'PAIR_SCORE_FAILURES_{label}.json').write_text(json.dumps(fails,indent=2)+'\n')
 (d/'PS02B_SCORE_SUMMARY.json').write_text(json.dumps({'kind':'PS02B_DEVELOPMENT_SCORES','formal_result':False,'model':name,'leaf_F_credit':{'win':1,'tie':0,'loss':-1},'formal_metric':'F_PS02B [-1,1]','panels':summ},indent=2)+'\n');print(json.dumps({'model':name,'panels':summ},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','esm2','prott5','esmc','score']);p.add_argument('--fasta',type=Path);p.add_argument('--checkpoint',type=Path);p.add_argument('--workdir',type=Path);a=p.parse_args()
 if a.command=='prepare':prepare()
 elif a.command=='score':score('esm2');score('prott5');score('esmc')
 else:
  if not all((a.fasta,a.checkpoint,a.workdir)):p.error('Forward requires --fasta, --checkpoint and --workdir')
  ROOT=a.workdir; ROOT.mkdir(parents=True,exist_ok=True); (ROOT/'endpoint_sequences.faa').write_bytes(a.fasta.read_bytes()); MODELS[a.command]=a.checkpoint
  encode(a.command)
