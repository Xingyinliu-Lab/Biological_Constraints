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
 out=[0]
 while out[-1]+W<L:
  nxt=min(out[-1]+stride,L-W)
  if nxt==out[-1]:break
  out.append(nxt)
 return out

def main():
 p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--weights',required=True);p.add_argument('--source-root',required=True);p.add_argument('--output',required=True);p.add_argument('--batch-size',type=int,default=4);p.add_argument('--device',default='cuda');a=p.parse_args();Path(a.output).parent.mkdir(parents=True,exist_ok=True);sys.path.insert(0,a.source_root);import fm
 model,alpha=fm.pretrained.rna_fm_t12(a.weights);model.eval().to(a.device);convert=alpha.get_batch_converter();rows=fasta(a.fasta);vals=[];ids=[];t=time.time();windowed=0
 def enc(batch):
  _,_,tok=convert(batch);h=model(tok.to(a.device),repr_layers=[12])['representations'][12]
  return [h[j,1:len(seq)+1].float() for j,(_,seq) in enumerate(batch)]
 with torch.inference_mode():
  for n,(ident,seq) in enumerate(rows):
   starts=layout(len(seq));windowed+=len(starts)>1;acc=None;counts=None
   chunks=[(ident,seq[s:s+1022]) for s in starts]
   for j in range(0,len(chunks),a.batch_size):
    for st,h in zip(starts[j:j+a.batch_size],enc(chunks[j:j+a.batch_size])):
     if acc is None:acc=torch.zeros((len(seq),h.shape[1]),device=a.device);counts=torch.zeros(len(seq),device=a.device)
     acc[st:st+h.shape[0]]+=h;counts[st:st+h.shape[0]]+=1
   if (counts==0).any():raise ValueError('uncovered nucleotide')
   vals.append((acc/counts[:,None]).mean(0).cpu().tolist());ids.append(ident)
   if n%100==0 or n+1==len(rows):print(json.dumps({'done':n+1,'total':len(rows),'windowed_nodes':windowed}),flush=True)
 z=np.asarray(vals,dtype=np.float32);order=np.argsort(np.array(ids));z=z[order];ids=np.array(ids)[order]
 if not np.isfinite(z).all():raise ValueError('nonfinite')
 np.savez(a.output,sequence_sha256=ids,vectors=z)
 Path(str(a.output)+'.metrics.json').write_text(json.dumps({'model_id':'RNA-FM','records':len(ids),'embedding_dim':int(z.shape[1]),'pooling':'overlapping 1022-nt coordinate windows; states averaged per nucleotide before equal-nucleotide whole-CDS mean; BOS/EOS excluded','window_nt':1022,'stride_nt':766,'windowed_nodes':windowed,'elapsed_seconds':time.time()-t,'weights_sha256':hashlib.sha256(Path(a.weights).read_bytes()).hexdigest(),'fasta_sha256':hashlib.sha256(Path(a.fasta).read_bytes()).hexdigest()},sort_keys=True)+'\n')
if __name__=='__main__':main()
