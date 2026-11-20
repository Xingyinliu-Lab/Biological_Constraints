from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np,torch

def fasta(p):
 r=[];k=None;s=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if k is not None:r.append((k,''.join(s)))
   k=x[1:].strip();s=[]
  else:s.append(x.strip().upper())
 if k is not None:r.append((k,''.join(s)))
 return r

def starts(L,W): return [0] if L<=W else [0,L-W]
def main():
 a=argparse.ArgumentParser();a.add_argument('--fasta',required=True);a.add_argument('--model-dir',required=True);a.add_argument('--output',required=True);a.add_argument('--device',default='cuda');z=a.parse_args();sys.path.insert(0,z.model_dir);from transformers import AutoModel,AutoTokenizer
 tok=AutoTokenizer.from_pretrained(z.model_dir,trust_remote_code=True);model=AutoModel.from_pretrained(z.model_dir,trust_remote_code=True,attn_implementation='eager').eval().to(z.device);nuc=set(tok.convert_tokens_to_ids(list('GAUCN')));out=[];ids=[];t=time.time();W=1022
 with torch.inference_mode():
  for ident,seq in fasta(z.fasta):
   acc=None;count=None
   for st in starts(len(seq),W):
    ch=seq[st:st+W];e=tok(ch,return_tensors='pt',padding=True,add_special_tokens=True);e={k:v.to(z.device) for k,v in e.items()};h=model(**e).last_hidden_state[0];m=torch.zeros_like(e['input_ids'][0],dtype=torch.bool)
    for q in nuc:m|=e['input_ids'][0].eq(q)
    hh=h[m].float()
    if hh.shape[0]!=len(ch):raise ValueError('nucleotide coordinate mismatch')
    if acc is None:acc=torch.zeros((len(seq),hh.shape[1]),device=z.device);count=torch.zeros(len(seq),device=z.device)
    acc[st:st+len(ch)]+=hh;count[st:st+len(ch)]+=1
   if (count==0).any():raise ValueError('uncovered nucleotide')
   out.append((acc/count[:,None]).mean(0).cpu().tolist());ids.append(ident)
 v=np.asarray(out,dtype=np.float32);np.savez(z.output,sequence_sha256=np.array(ids),vectors=v)
 Path(str(z.output)+'.metrics.json').write_text(json.dumps({'model_id':'ERNIE-RNA','records':len(ids),'window_nt':W,'layout':'two coordinate-anchored overlapping windows for boundary sequences; overlap states averaged per nucleotide before whole-object mean','elapsed_seconds':time.time()-t,'weights_sha256':hashlib.sha256(Path(z.model_dir,'pytorch_model.bin').read_bytes()).hexdigest()},sort_keys=True)+'\n')
if __name__=='__main__':main()