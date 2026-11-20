#!/usr/bin/env python3
"""RNA/PS09-compatible ESM3 sequence-only residue endpoint extractor."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,time
from collections import defaultdict
from pathlib import Path
import numpy as np, torch
from safetensors.numpy import save_file
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein
WINDOW=1023; STRIDE=767; INTERFACE='ESM3_SEQ_MODEL_INTERFACE_v1.0.0_FROZEN_20261010'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def fasta(p):
 d={}; k=None; b=[]
 for q in open(p):
  q=q.strip()
  if not q: continue
  if q[0]=='>':
   if k: d[k]=''.join(b)
   k=q[1:].split()[0]; b=[]
  else:b.append(q.upper())
 if k:d[k]=''.join(b)
 return d
def layout(n):
 if n<=WINDOW:return [0]
 a=list(range(0,n-WINDOW+1,STRIDE))
 if a[-1]!=n-WINDOW:a.append(n-WINDOW)
 return a
def eid(r):return f"{r['block_id']}|{r['state']}|{r['uniprot_position_1based']}"
def main():
 a=argparse.ArgumentParser();a.add_argument('--membership',type=Path,required=True);a.add_argument('--fasta',type=Path,required=True);a.add_argument('--output-dir',type=Path,required=True);x=a.parse_args();x.output_dir.mkdir(parents=True,exist_ok=False)
 with gzip.open(x.membership,'rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
 seq=fasta(x.fasta); positions=defaultdict(list); seen=set()
 for r in rows:
  k=eid(r); pos=int(r['uniprot_position_1based'])-1; p=r['uniprot_accession']
  if k in seen or p not in seq or not 0<=pos<len(seq[p]) or seq[p][pos]!=r['amino_acid']:raise RuntimeError('frozen endpoint closure')
  seen.add(k);positions[p].append((pos,k))
 dev=torch.device('cuda:0');torch.cuda.set_device(dev);torch.cuda.reset_peak_memory_stats(dev);m=ESM3.from_pretrained('esm3_sm_open_v1',device=dev).eval();out={};idx=[];tile_count=0;t0=time.time()
 with torch.inference_mode(),torch.autocast('cuda',dtype=torch.bfloat16):
  for p,s in seq.items():
   if p not in positions:continue
   acc=None;count=None
   for left in layout(len(s)):
    right=min(len(s),left+WINDOW); t=m.encode(ESMProtein(sequence=s[left:right])).sequence.unsqueeze(0).to(dev);h=m.forward(sequence_tokens=t).embeddings[0,1:-1].float().cpu()
    if acc is None:acc=torch.zeros((len(s),h.shape[1]));count=torch.zeros(len(s))
    acc[left:right]+=h;count[left:right]+=1;tile_count+=1
   if (count==0).any():raise RuntimeError('coverage')
   h=acc/count[:,None]
   for pos,k in positions[p]:
    v=h[pos].numpy().astype(np.float32);out[k]=v;idx.append((k,p,pos,k,len(v),hashlib.sha256(v.tobytes()).hexdigest()))
 if set(out)!=seen:raise RuntimeError('endpoint output closure')
 dest=x.output_dir/'endpoint_embeddings.safetensors';save_file(out,str(dest))
 ip=x.output_dir/'embedding_index.tsv'
 with ip.open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(('endpoint_id','uniprot_accession','position_0based','embedding_key','embedding_dim','embedding_sha256'));w.writerows(idx)
 rep={'status':'PASS','model_interface_release_id':INTERFACE,'track_policy':'sequence_only; non-sequence tracks omitted and default masked','sequence_window_residues':WINDOW,'sequence_stride_residues':STRIDE,'compute_dtype':'bfloat16','stored_embedding_dtype':'float32','endpoint_count':len(idx),'protein_count':len(positions),'tile_count':tile_count,'embedding_dim':len(idx) and idx[0][4],'membership_sha256':sha(x.membership),'fasta_sha256':sha(x.fasta),'embedding_sha256':sha(dest),'embedding_index_sha256':sha(ip),'elapsed_seconds':time.time()-t0,'peak_gpu_mib':torch.cuda.max_memory_allocated(dev)/2**20}
 (x.output_dir/'FORWARD_REPORT.json').write_text(json.dumps(rep,sort_keys=True,indent=2)+'\n')
if __name__=='__main__':main()
