from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np,torch
from gb_rna import RNABertModel,RNABertTokenizer,get_vocab_filepath,validate_sequences

def fasta(p):
 r=[];i=None;s=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if i is not None:r.append((i,''.join(s)))
   i=x[1:].strip();s=[]
  else:s.append(x.strip().upper())
 if i is not None:r.append((i,''.join(s)))
 for ident,seq in r:
  if hashlib.sha256(seq.encode()).hexdigest()!=ident: raise ValueError('sequence SHA256 identity mismatch')
 if len({x[0] for x in r})!=len(r): raise ValueError('duplicate sequence identity')
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
 p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--model-dir',required=True);p.add_argument('--output',required=True);p.add_argument('--batch-size',type=int,default=2);p.add_argument('--device',default='cuda');a=p.parse_args();rows=fasta(a.fasta);tok=RNABertTokenizer(get_vocab_filepath(),version='v2');model=RNABertModel.from_pretrained(a.model_dir).to(a.device).to(torch.bfloat16).eval();valid=torch.tensor([5,6,7,8,9,10],device=a.device);vals={};t=time.time();windowed=0
 def enc(chunks):
  q=tok(validate_sequences(chunks),max_length=1024,padding=True,truncation=False,return_tensors='pt');q={k:v.to(a.device) for k,v in q.items()};h=model(input_ids=q['input_ids'],attention_mask=q['attention_mask']).last_hidden_state.float();m=(q['input_ids'].unsqueeze(-1)==valid).any(-1)
  if (m.sum(1)==0).any():raise ValueError('no valid nucleotide token')
  return [h[k][m[k]] for k in range(len(chunks))]
 with torch.inference_mode():
  for n,(ident,seq) in enumerate(rows):
   starts=layout(len(seq));windowed+=len(starts)>1;acc=None;counts=None;chunks=[seq[s:s+1022] for s in starts]
   for j in range(0,len(chunks),a.batch_size):
    for st,h in zip(starts[j:j+a.batch_size],enc(chunks[j:j+a.batch_size])):
     if h.shape[0]!=len(seq[st:st+1022]):raise ValueError('nucleotide coordinate mismatch')
     if acc is None:acc=torch.zeros((len(seq),h.shape[1]),device=a.device);counts=torch.zeros(len(seq),device=a.device)
     acc[st:st+h.shape[0]]+=h;counts[st:st+h.shape[0]]+=1
   if (counts==0).any():raise ValueError('uncovered nucleotide')
   vals[ident]=(acc/counts[:,None]).float().cpu()
   if n%50==0 or n+1==len(rows):print(json.dumps({'done':n+1,'total':len(rows),'windowed_nodes':windowed}),flush=True)
 z=torch.cat([v.reshape(-1) for v in vals.values()])
 if not torch.isfinite(z).all():raise ValueError('nonfinite')
 torch.save(vals,a.output)
 Path(str(a.output)+'.metrics.json').write_text(json.dumps({'model_id':'AIDO/GB.RNA-650M','records':len(vals),'representation':'per-nucleotide hidden state vectors','pooling':'overlapping 1022-nt coordinate windows; states averaged per nucleotide before equal-nucleotide whole-CDS mean; special tokens excluded','window_nt':1022,'stride_nt':766,'windowed_nodes':windowed,'elapsed_seconds':time.time()-t,'weights_sha256':hashlib.sha256((Path(a.model_dir)/'model.safetensors').read_bytes()).hexdigest()},sort_keys=True)+'\n')
if __name__=='__main__':main()