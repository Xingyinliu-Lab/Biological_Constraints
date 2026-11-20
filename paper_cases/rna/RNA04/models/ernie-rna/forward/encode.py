from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np,torch

def fasta(p):
 r=[];i=None;s=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if i is not None:r.append((i,''.join(s)))
   i=x[1:].strip();s=[]
  else:s.append(x.strip().upper())
 if i is not None:r.append((i,''.join(s)))
 return r

def layout(L,W=1022,stride=766):
 if L<=W:return [0]
 o=[0]
 while o[-1]+W<L:
  n=min(o[-1]+stride,L-W)
  if n==o[-1]:break
  o.append(n)
 return o

def main():
 p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--model-dir',required=True);p.add_argument('--output',required=True);p.add_argument('--device',default='cuda');a=p.parse_args();Path(a.output).parent.mkdir(parents=True,exist_ok=True);sys.path.insert(0,a.model_dir);from transformers import AutoModel,AutoTokenizer
 tok=AutoTokenizer.from_pretrained(a.model_dir,trust_remote_code=True);model=AutoModel.from_pretrained(a.model_dir,trust_remote_code=True,attn_implementation='eager').eval().to(a.device);nuc=set(tok.convert_tokens_to_ids(list('GAUCN')));rows=fasta(a.fasta);vals=[];ids=[];t=time.time();windowed=0
 with torch.inference_mode():
  for n,(ident,seq) in enumerate(rows):
   starts=layout(len(seq));windowed+=len(starts)>1;acc=None;counts=None
   for st in starts:
    ch=seq[st:st+1022];e=tok(ch,return_tensors='pt',padding=True,add_special_tokens=True);e={k:v.to(a.device) for k,v in e.items()};h=model(**e).last_hidden_state[0];m=torch.zeros_like(e['input_ids'][0],dtype=torch.bool)
    for q in nuc:m|=e['input_ids'][0].eq(q)
    hh=h[m].float()
    if hh.shape[0]!=len(ch):raise ValueError('nucleotide coordinate mismatch')
    if acc is None:acc=torch.zeros((len(seq),hh.shape[1]),device=a.device);counts=torch.zeros(len(seq),device=a.device)
    acc[st:st+len(ch)]+=hh;counts[st:st+len(ch)]+=1
   if (counts==0).any():raise ValueError('uncovered nucleotide')
   vals.append((acc/counts[:,None]).mean(0).cpu().tolist());ids.append(ident)
   if n%50==0 or n+1==len(rows):print(json.dumps({'done':n+1,'total':len(rows),'windowed_nodes':windowed}),flush=True)
 z=np.asarray(vals,dtype=np.float32)
 if not np.isfinite(z).all():raise ValueError('nonfinite')
 np.savez(a.output,sequence_sha256=np.array(ids),vectors=z)
 Path(str(a.output)+'.metrics.json').write_text(json.dumps({'model_id':'ERNIE-RNA','records':len(ids),'embedding_dim':int(z.shape[1]),'pooling':'overlapping 1022-nt coordinate windows; states averaged per nucleotide before equal-nucleotide whole-CDS mean; CLS/EOS/pad excluded','window_nt':1022,'stride_nt':766,'windowed_nodes':windowed,'attention':'eager full-context ERNIE-RNA compatibility checkpoint','elapsed_seconds':time.time()-t,'weights_sha256':hashlib.sha256(Path(a.model_dir,'pytorch_model.bin').read_bytes()).hexdigest()},sort_keys=True)+'\n')
if __name__=='__main__':main()
