#!/usr/bin/env python3
"""Convert any complete PS02B representation set into formal leaf scores."""
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path
import numpy as np
from ps02_safe_cosine_total_float32 import safe_cosine_total_float32

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def rows(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta_ids(p):
 ids=[]
 with open(p,'rt') as f:
  for line in f:
   if line.startswith('>'):ids.append(line[1:].split()[0])
 if not ids or len(ids)!=len(set(ids)):raise ValueError('invalid or duplicate FASTA IDs')
 return set(ids)
def main():
 p=argparse.ArgumentParser()
 for x in ('relations','embeddings','output','model-id'):p.add_argument('--'+x,required=True)
 p.add_argument('--embedding-index',help='Required only for Safetensors input')
 p.add_argument('--objects-fasta',required=True,help='Released Forward FASTA used to validate exact object coverage')
 p.add_argument('--binding',help='Optional run manifest written for the supplied model')
 p.add_argument('--checkpoint-id',help='Optional checkpoint/revision identifier recorded in --binding')
 a=p.parse_args(); rel=rows(a.relations)
 used={x for r in rel for x in (r['anchor'],r['positive'],r['negative'])}
 if str(a.embeddings).lower().endswith('.npz'):
  with np.load(a.embeddings,allow_pickle=False) as archive:
   if set(archive.files)!={'object_ids','vectors'}:raise ValueError('NPZ must contain exactly object_ids and vectors')
   ids=[str(x) for x in archive['object_ids']];matrix=np.asarray(archive['vectors'],dtype=np.float32)
  if matrix.ndim!=2 or matrix.shape[1]<=0 or matrix.shape[0]!=len(ids) or len(ids)!=len(set(ids)):raise ValueError('invalid NPZ shape or duplicate object ID')
  if set(ids)!=fasta_ids(a.objects_fasta):raise ValueError('representation IDs do not exactly match released Forward FASTA')
  if not used<=set(ids):raise ValueError('formal scored endpoint absent from representations')
  if not np.isfinite(matrix).all():raise ValueError('non-finite representation')
  emb={x:matrix[i] for i,x in enumerate(ids)}
 else:
  if not a.embedding_index:raise ValueError('--embedding-index is required for Safetensors input')
  from safetensors import safe_open
  index_rows=rows(a.embedding_index)
  ids=[r['domain_id'] for r in index_rows];keys=[r['embedding_key'] for r in index_rows]
  if len(ids)!=len(set(ids)) or len(keys)!=len(set(keys)):raise ValueError('duplicate domain_id or embedding_key')
  idx=dict(zip(ids,keys))
  if set(idx)!=fasta_ids(a.objects_fasta):raise ValueError('embedding index does not exactly match released Forward FASTA')
  if not used <= set(idx):raise ValueError('formal scored endpoint absent from embedding index')
  with safe_open(a.embeddings,framework='np') as fh:
   if set(fh.keys())!=set(keys):raise ValueError('Safetensors keys do not match embedding index')
   emb={x:np.asarray(fh.get_tensor(idx[x]),dtype=np.float32) for x in used}
  if any(v.ndim!=1 or not np.isfinite(v).all() for v in emb.values()):raise ValueError('invalid Safetensors representations')
  dims={v.shape for v in emb.values()}
  if len(dims)!=1 or next(iter(dims))[0]<=0:raise ValueError('inconsistent or zero-dimensional Safetensors representations')
 out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
 with gzip.open(out,'wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['bundle_id','score_anchor_positive','score_anchor_negative'],delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rel:w.writerow({'bundle_id':r['bundle_id'],'score_anchor_positive':repr(safe_cosine_total_float32(emb[r['anchor']],emb[r['positive']])),'score_anchor_negative':repr(safe_cosine_total_float32(emb[r['anchor']],emb[r['negative']]))})
 if a.binding:
  producer=Path(__file__).resolve()
  manifest={'model_id':a.model_id,'leaf_score_producer_sha256':sha(producer),'safe_cosine_operator_id':'PS02_SAFE_COSINE_TOTAL_FLOAT32','safe_cosine_implementation_sha256':sha(Path(__file__).with_name('ps02_safe_cosine_total_float32.py')),'leaf_score_ledger_sha256':sha(out)}
  if a.checkpoint_id:manifest['checkpoint_id']=a.checkpoint_id
  binding=Path(a.binding);binding.parent.mkdir(parents=True,exist_ok=True);binding.write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
