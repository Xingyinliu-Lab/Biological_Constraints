#!/usr/bin/env python3
"""PS04 S/L model-independent endpoint encoding and schema-bound scoring."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math,re,sys,time
from pathlib import Path
import numpy as np
import torch
from safetensors.torch import save_file
from safetensors import safe_open
BASE=Path.cwd(); ROOT=None; REL=None; FASTA=None
PANELS={'S_2K':REL/'S_2K','L_4K':REL/'L_4K'}
MODELS={'esm2':None,'prott5':None,'esmc':None}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for q in iter(lambda:f.read(1<<20),b''):h.update(q)
 return h.hexdigest()
def gzrows(p):
 with gzip.open(p,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def objects():
 return {r['domain_id']:r for r in gzrows(PANELS['L_4K']/'OBJECT_LEDGER.tsv.gz')}
def seqs(want):
 out={}; name=None; chunks=[]
 with FASTA.open() as f:
  for line in f:
   if line.startswith('>'):
    if name in want:out[name]=''.join(chunks)
    name=line[1:].split('|')[-1].split('/')[0];chunks=[]
   else:chunks.append(line.strip())
  if name in want:out[name]=''.join(chunks)
 return out
def read_fasta(p):
 out={};k=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):k=line[1:];out[k]=''
  else:out[k]+=line.strip()
 return out
def batches(s,n,total):
 batch=[];size=0
 for i in sorted(s,key=lambda x:(len(s[x]),x)):
  z=len(s[i])+2
  if batch and (len(batch)>=n or size+z>total):yield batch;batch=[];size=0
  batch.append(i);size+=z
 if batch:yield batch
def prepare():
 rel={k:gzrows(v/'LEAF_RELATION_LEDGER.tsv.gz') for k,v in PANELS.items()}
 ids={x for z in rel.values() for r in z for x in (r['anchor_id'],r['positive_id'],r['negative_id'])}
 obj=objects(); s=seqs(ids); bad=[]
 for i in sorted(ids):
  if i not in obj:bad.append([i,'OBJECT_MISSING'])
  elif i not in s:bad.append([i,'SEQUENCE_MISSING'])
  elif hashlib.sha256(s[i].encode()).hexdigest()!=obj[i]['sequence_sha256']:bad.append([i,'SEQUENCE_HASH_MISMATCH'])
 ROOT.mkdir(parents=True,exist_ok=True)
 (ROOT/'endpoint_sequences.faa').write_text(''.join(f'>{i}\n{s[i]}\n' for i in sorted(s)))
 rep={'kind':'PS04_ENCODING_INPUT','status':'PASS' if not bad else 'FAIL','endpoint_count':len(ids),'sequence_found':len(s),'failure_count':len(bad),'failures':bad[:100],'source_fasta_sha256':sha(FASTA),'panels':{k:{'leaf_count':len(v),'leaf_sha256':sha(PANELS[k]/'LEAF_RELATION_LEDGER.tsv.gz')} for k,v in rel.items()}}
 (ROOT/'PREPARE_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));
 if bad:raise SystemExit('PREPARE_BLOCKED')
def encode(name):
 s=read_fasta(ROOT/'endpoint_sequences.faa');out=ROOT/name;out.mkdir(exist_ok=True);dest=out/'embeddings.safetensors'
 if dest.exists():raise SystemExit('EXISTING_EMBEDDINGS')
 torch.cuda.set_device(0);torch.cuda.reset_peak_memory_stats();torch.manual_seed(20261001)
 if name=='esm2':
  torch.serialization.add_safe_globals([argparse.Namespace])
  sys.path.insert(0,str(BASE/'m3_smoke/python_deps'));import esm
  m,a=esm.pretrained.load_model_and_alphabet_local(str(MODELS[name]));m=m.eval().half().cuda();conv=a.get_batch_converter();meta={'model':'esm2_t33_650M_UR50D','layer':33,'pooling':'residue_mean'};bi=batches(s,16,4096)
  def f(ids):
   _,_,tok=conv([(i,s[i]) for i in ids])
   with torch.inference_mode():h=m(tok.cuda(),repr_layers=[33],return_contacts=False)['representations'][33]
   return [h[j,1:len(s[i])+1].float().mean(0).cpu().contiguous() for j,i in enumerate(ids)]
 elif name=='prott5':
  sys.path.insert(0,str(BASE/'m3_encoder_audit/python_deps_prott5_444'));from transformers import T5EncoderModel,T5Tokenizer
  tok=T5Tokenizer.from_pretrained(str(MODELS[name]),do_lower_case=False,legacy=True,local_files_only=True);m=T5EncoderModel.from_pretrained(str(MODELS[name]),torch_dtype=torch.float16,local_files_only=True).cuda().eval();meta={'model':'Rostlab/prot_t5_xl_half_uniref50-enc','pooling':'residue_mean','UO_mapping':'X'};bi=batches(s,8,2048)
  def f(ids):
   enc=tok([' '.join(re.sub(r'[UZOB]','X',s[i])) for i in ids],add_special_tokens=True,padding=True,return_tensors='pt')
   with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):h=m(input_ids=enc['input_ids'].cuda(),attention_mask=enc['attention_mask'].cuda()).last_hidden_state
   return [h[j,:len(s[i])].float().mean(0).cpu().contiguous() for j,i in enumerate(ids)]
 else:
  from esm.models.esmc import EsmcForMaskedLM,EsmcTokenizer
  tok=EsmcTokenizer();m=EsmcForMaskedLM.from_pretrained(str(MODELS[name]),local_files_only=True).cuda().eval();meta={'model':'biohub/ESMC-600M','pooling':'residue_mean'};bi=batches(s,16,4096)
  def f(ids):
   enc={k:v.cuda() for k,v in tok([s[i] for i in ids],return_tensors='pt',padding=True).items()}
   with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):h=m(**enc).last_hidden_state
   return [h[j,1:len(s[i])+1].float().mean(0).cpu().contiguous() for j,i in enumerate(ids)]
 vec={};idx=[];fails=[];start=time.time()
 for nb,ids in enumerate(bi,1):
  try:r=f(ids)
  except Exception as e:
   r=[]
   for i in ids:
    try:r.append((i,f([i])[0]))
    except Exception as e2:fails.append([i,type(e2).__name__,str(e2)[:160]])
   for i,x in r:
    if torch.isfinite(x).all():
     k=f'embedding_{len(vec)+1:06d}';vec[k]=x;idx.append([i,k,len(s[i]),x.numel(),hashlib.sha256(x.numpy().tobytes()).hexdigest()])
   continue
  for i,x in zip(ids,r):
   if torch.isfinite(x).all():
    k=f'embedding_{len(vec)+1:06d}';vec[k]=x;idx.append([i,k,len(s[i]),x.numel(),hashlib.sha256(x.numpy().tobytes()).hexdigest()])
   else:fails.append([i,'NONFINITE',''])
  if nb%50==0:print(json.dumps({'model':name,'batches':nb,'encoded':len(idx),'failed':len(fails)}),flush=True)
 save_file(vec,str(dest),metadata={k:str(v) for k,v in meta.items()})
 with (out/'embedding_index.tsv').open('w',newline='') as q:
  w=csv.writer(q,delimiter='\t',lineterminator='\n');w.writerow(['domain_id','embedding_key','sequence_length','embedding_dim','embedding_sha256']);w.writerows(idx)
 with (out/'failure_ledger.tsv').open('w',newline='') as q:
  w=csv.writer(q,delimiter='\t',lineterminator='\n');w.writerow(['domain_id','failure_type','detail']);w.writerows(fails)
 rep={'kind':'PS04_MODEL_FORWARD','status':'PASS' if len(idx)==len(s) and not fails else 'PARTIAL','model':name,'endpoint_count':len(s),'encoded_count':len(idx),'failure_count':len(fails),'checkpoint_sha256':sha(MODELS[name] if MODELS[name].is_file() else MODELS[name]/('pytorch_model.bin' if name=='prott5' else 'model.safetensors')),'embedding_sha256':sha(dest),'embedding_dim':idx[0][3] if idx else None,'seconds':time.time()-start,'peak_gpu_mib':torch.cuda.max_memory_allocated()/2**20}
 (out/'FORWARD_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep),flush=True)
def score(name):
 idx={r['domain_id']:r['embedding_key'] for r in csv.DictReader((ROOT/name/'embedding_index.tsv').open(),delimiter='\t')}
 with safe_open(str(ROOT/name/'embeddings.safetensors'),framework='np') as f:
  ids=sorted(idx);E=np.stack([f.get_tensor(idx[i]).astype(np.float64) for i in ids]);loc={x:n for n,x in enumerate(ids)}
 E/=np.maximum(np.linalg.norm(E,axis=1,keepdims=True),1e-300)
 for panel,path in PANELS.items():
  out=ROOT/name/f'LEAF_SCORES_{panel}.tsv'; n=0
  with out.open('w',newline='') as q:
   fields=['negative_endpoint_id','negative_hsf_id','negative_topology_id','positive_endpoint_id','positive_hsf_id','anchor_id','anchor_parent_id','anchor_source_unit_id','anchor_hsf_id','anchor_topology_id','anchor_architecture_id','positive_score','negative_score'];w=csv.DictWriter(q,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
   rows=[]
   for r in gzrows(path/'LEAF_RELATION_LEDGER.tsv.gz'):
    rows.append(r)
    if len(rows)>=50000:
     emit(rows,w,loc,E);n+=len(rows);rows=[]
   if rows:emit(rows,w,loc,E);n+=len(rows)
  print(json.dumps({'model':name,'panel':panel,'leaf_scores':n}),flush=True)
def emit(rows,w,loc,E):
 ai=np.fromiter((loc[r['anchor_id']] for r in rows),int);pi=np.fromiter((loc[r['positive_id']] for r in rows),int);ni=np.fromiter((loc[r['negative_id']] for r in rows),int)
 ps=np.einsum('ij,ij->i',E[ai],E[pi]);ns=np.einsum('ij,ij->i',E[ai],E[ni])
 for r,x,y in zip(rows,ps,ns):w.writerow({'negative_endpoint_id':r['negative_id'],'negative_hsf_id':r['negative_hsf_id'],'negative_topology_id':r['negative_topology_id'],'positive_endpoint_id':r['positive_id'],'positive_hsf_id':r['positive_hsf_id'],'anchor_id':r['anchor_id'],'anchor_parent_id':r['anchor_parent_proxy'],'anchor_source_unit_id':r['anchor_source_unit_id'],'anchor_hsf_id':r['anchor_hsf_id'],'anchor_topology_id':r['anchor_topology_id'],'anchor_architecture_id':r['anchor_architecture_id'],'positive_score':repr(float(x)),'negative_score':repr(float(y))})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('cmd',choices=['prepare','esm2','prott5','esmc','score']);p.add_argument('--fasta',type=Path);p.add_argument('--checkpoint',type=Path);p.add_argument('--workdir',type=Path);a=p.parse_args()
 if any((a.fasta,a.checkpoint,a.workdir)):
  if a.cmd not in MODELS or not all((a.fasta,a.checkpoint,a.workdir)):p.error('current Forward requires a model command plus --fasta, --checkpoint and --workdir')
  ROOT=a.workdir;ROOT.mkdir(parents=True,exist_ok=True);(ROOT/'endpoint_sequences.faa').write_bytes(a.fasta.read_bytes());MODELS[a.cmd]=a.checkpoint
 if a.cmd=='prepare':prepare()
 elif a.cmd=='score':
  for n in ('esm2','prott5','esmc'):score(n)
 else:encode(a.cmd)
