#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, time
from pathlib import Path
import torch
from safetensors.torch import save_file
from sequence_models.pretrained import load_model_and_alphabet
ROOT=None
INPUT=None
RUNTIME=None
CHECKPOINT=None
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(1048576),b''): h.update(c)
 return h.hexdigest()
def read_fasta(p):
 out={}; k=None; a=[]
 for line in p.read_text().splitlines():
  if line.startswith('>'):
   if k is not None: out[k]=''.join(a)
   k=line[1:].split()[0]; a=[]
  else: a.append(line.strip())
 if k is not None: out[k]=''.join(a)
 return out
def batches(keys, seqs, n=16, cap=4096):
 cur=[]; total=0
 for k in keys:
  l=len(seqs[k])
  if cur and (len(cur)>=n or total+l>cap):
   yield cur; cur=[]; total=0
  cur.append(k); total+=l
 if cur: yield cur
def main():
 global ROOT, INPUT, RUNTIME, CHECKPOINT
 parser=argparse.ArgumentParser()
 parser.add_argument('--fasta',type=Path,required=True)
 parser.add_argument('--checkpoint',type=Path,required=True)
 parser.add_argument('--workdir',type=Path,required=True)
 args=parser.parse_args()
 if any((args.fasta,args.checkpoint,args.workdir)):
  if not all((args.fasta,args.checkpoint,args.workdir)):parser.error('current Forward requires --fasta, --checkpoint and --workdir')
  INPUT=args.fasta;CHECKPOINT=args.checkpoint;ROOT=args.workdir;RUNTIME=CHECKPOINT.parents[2]
 os.environ.setdefault('TORCH_HOME',str(CHECKPOINT.parents[2]))
 out=ROOT/'carp'; out.mkdir(parents=True,exist_ok=True)
 dest=out/'embeddings.safetensors'
 if dest.exists(): raise SystemExit('EXISTING_OUTPUT_PRESERVED')
 seqs=read_fasta(INPUT); keys=list(seqs)
 if len(keys)!=13720 or len(keys)!=len(set(keys)): raise ValueError('INPUT_POPULATION_MISMATCH')
 torch.cuda.set_device(0); torch.cuda.reset_peak_memory_stats(); torch.manual_seed(20261001)
 model, alphabet=load_model_and_alphabet('carp_640M'); model=model.eval().cuda()
 layer=len(model.model.embedder.layers)
 def forward(ids):
  tokens=alphabet([(seqs[x],) for x in ids])[0].cuda()
  with torch.inference_mode(), torch.autocast('cuda',dtype=torch.float16):
   h=model(tokens,repr_layers=[-1],logits=False)['representations'][layer]
  return [h[j,:len(seqs[i])].float().mean(0).cpu().contiguous() for j,i in enumerate(ids)]
 vectors={}; idx=[]; failures=[]; started=time.time()
 for num,ids in enumerate(batches(keys,seqs),1):
  try: result=forward(ids)
  except Exception as e:
   if len(ids)>1:
    for ident in ids:
     try: result1=forward([ident])[0]
     except Exception as x: failures.append((ident,type(x).__name__,str(x)[:180])); continue
     if not torch.isfinite(result1).all(): failures.append((ident,'NONFINITE_EMBEDDING','')); continue
     key=f'embedding_{len(vectors)+1:06d}'; vectors[key]=result1; idx.append((ident,key,len(seqs[ident]),result1.numel(),hashlib.sha256(result1.numpy().tobytes()).hexdigest()))
    continue
   failures.append((ids[0],type(e).__name__,str(e)[:180])); continue
  for ident,vec in zip(ids,result):
   if not torch.isfinite(vec).all(): failures.append((ident,'NONFINITE_EMBEDDING','')); continue
   key=f'embedding_{len(vectors)+1:06d}'; vectors[key]=vec; idx.append((ident,key,len(seqs[ident]),vec.numel(),hashlib.sha256(vec.numpy().tobytes()).hexdigest()))
  if num%50==0: print(json.dumps({'batches':num,'encoded':len(idx),'failed':len(failures)}),flush=True)
 if vectors: save_file(vectors,str(dest),metadata={'model':'CARP-640M','architecture':'ByteNet_CNN','pooling':'residue_mean','representation_layer':str(layer)})
 with (out/'embedding_index.tsv').open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(('occurrence_id','embedding_key','sequence_length','embedding_dim','embedding_sha256'));w.writerows(idx)
 with (out/'failure_ledger.tsv').open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(('occurrence_id','failure_type','detail'));w.writerows(failures)
 report={'status':'PASS' if len(idx)==len(seqs) and not failures else 'PARTIAL','kind':'PS02A_DEVELOPMENT_MODEL_FORWARD','formal_result':False,'formal_model_interface_release':False,'model':'CARP-640M','architecture':'ByteNet_CNN','pooling':'RESIDUE_MEAN','representation_layer':layer,'endpoint_count':len(seqs),'encoded_count':len(idx),'failure_count':len(failures),'checkpoint_sha256':sha(CHECKPOINT),'embedding_file_sha256':sha(dest) if dest.exists() else None,'embedding_dim':idx[0][3] if idx else None,'seconds':time.time()-started,'peak_gpu_memory_mib':torch.cuda.max_memory_allocated()/2**20}
 (out/'FORWARD_REPORT.json').write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report),flush=True)
if __name__=='__main__': main()
