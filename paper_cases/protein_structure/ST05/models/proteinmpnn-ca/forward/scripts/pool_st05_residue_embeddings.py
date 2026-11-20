#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def mean_embedding(p):
    with np.load(p,allow_pickle=False) as z:
        x=np.asarray(z['residue_embedding'],dtype=np.float64)
    if x.ndim!=2 or x.shape[0]==0 or not np.isfinite(x).all(): raise ValueError(f'invalid residue embedding: {p}')
    return x.mean(axis=0,dtype=np.float64).astype(np.float32),x.shape[0]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--relations',type=Path,required=True); ap.add_argument('--residue-dir',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--model-id',required=True); a=ap.parse_args()
    with a.relations.open(encoding='utf-8',newline='') as h: rows=list(csv.DictReader(h,delimiter='\t'))
    bases=sorted({r['base_object_id'] for r in rows}); relids=[r['relation_id'] for r in rows]
    if len(relids)!=len(set(relids)): raise ValueError('duplicate relation')
    bvec=[]; vvec=[]; counts={}; records=[]
    for oid in bases:
        p=a.residue_dir/f'{oid}__BASE.npz'; v,n=mean_embedding(p); bvec.append(v); counts[oid]=n; records.append({'id':oid,'role':'base','residues':n,'sha256':sha(p)})
    for r in rows:
        oid=r['variant_object_id']; p=a.residue_dir/f'{oid}.npz'; v,n=mean_embedding(p); vvec.append(v); records.append({'id':r['relation_id'],'variant_object_id':oid,'role':'variant','residues':n,'sha256':sha(p)})
        if n!=counts[r['base_object_id']]: raise ValueError(f'residue-count closure failure: {oid}')
    B=np.asarray(bvec,dtype=np.float32); V=np.asarray(vvec,dtype=np.float32)
    if B.ndim!=2 or V.ndim!=2 or B.shape[1]!=V.shape[1]: raise ValueError('dimension closure failure')
    a.output.parent.mkdir(parents=True,exist_ok=True); np.savez_compressed(a.output,base_object_ids=np.asarray(bases),base_embeddings=B,relation_ids=np.asarray(relids),variant_embeddings=V)
    rep={'status':'PASS','model_id':a.model_id,'base_objects':len(bases),'relations':len(rows),'dimension':int(B.shape[1]),'pooling':'float64 residue mean then float32','embedding_sha256':sha(a.output),'records':records}
    a.output.with_name('FORWARD_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps({k:rep[k] for k in ('status','model_id','base_objects','relations','dimension','embedding_sha256')},sort_keys=True))
if __name__=='__main__': main()
