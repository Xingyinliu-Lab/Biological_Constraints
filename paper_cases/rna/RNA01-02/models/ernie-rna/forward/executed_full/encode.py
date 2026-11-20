from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np, torch

def fasta(path):
 out=[]; ident=None; seq=[]
 for line in Path(path).read_text().splitlines():
  if line.startswith('>'):
   if ident is not None: out.append((ident,''.join(seq)))
   ident=line[1:].strip(); seq=[]
  else: seq.append(line.strip().upper())
 if ident is not None: out.append((ident,''.join(seq)))
 return out

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--fasta',required=True); ap.add_argument('--model-dir',required=True); ap.add_argument('--output',required=True); ap.add_argument('--batch-size',type=int,default=1); ap.add_argument('--limit',type=int); ap.add_argument('--device',default='cuda'); a=ap.parse_args()
 sys.path.insert(0,a.model_dir)
 from transformers import AutoModel, AutoTokenizer
 tok=AutoTokenizer.from_pretrained(a.model_dir,trust_remote_code=True)
 model=AutoModel.from_pretrained(a.model_dir,trust_remote_code=True,attn_implementation='eager').eval().to(a.device)
 rows=fasta(a.fasta); rows=rows[:a.limit] if a.limit else rows; vals=[]; ids=[]; start=time.time(); window=1022; split_nodes=0
 nuc_ids=set(tok.convert_tokens_to_ids(list('GAUCN')))
 with torch.inference_mode():
  for i,(ident,seq) in enumerate(rows):
   chunks=[seq[k:k+window] for k in range(0,len(seq),window)]
   if len(chunks)>1: split_nodes+=1
   states=[]
   for chunk in chunks:
    enc=tok(chunk,return_tensors='pt',padding=True,add_special_tokens=True); enc={k:v.to(a.device) for k,v in enc.items()}
    h=model(**enc).last_hidden_state[0]
    mask=torch.zeros_like(enc['input_ids'][0],dtype=torch.bool)
    for q in nuc_ids: mask |= enc['input_ids'][0].eq(q)
    states.append(h[mask])
   z=torch.cat(states,dim=0).mean(dim=0)
   vals.append(z.float().cpu().tolist()); ids.append(ident)
   if i%50==0 or i+1==len(rows): print(json.dumps({'done':i+1,'total':len(rows),'windowed_boundary_nodes':split_nodes}),flush=True)
 z=np.asarray(vals,dtype=np.float32)
 if not np.isfinite(z).all(): raise ValueError('nonfinite')
 np.savez(a.output,sequence_sha256=np.array(ids),vectors=z)
 Path(str(a.output)+'.metrics.json').write_text(json.dumps({'model_id':'ERNIE-RNA','records':len(ids),'embedding_dim':int(z.shape[1]),'pooling':'mean valid nucleotide token states across whole object; CLS/EOS/pad excluded','windowed_boundary_nodes':split_nodes,'window_nt':window,'attention':'eager full-context ERNIE-RNA compatibility checkpoint','weights_sha256':hashlib.sha256(Path(a.model_dir,'pytorch_model.bin').read_bytes()).hexdigest(),'elapsed_seconds':time.time()-start},sort_keys=True)+'\n')
if __name__=='__main__': main()
