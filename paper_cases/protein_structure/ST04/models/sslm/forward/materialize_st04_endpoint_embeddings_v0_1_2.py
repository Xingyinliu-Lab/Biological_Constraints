#!/usr/bin/env python3
import argparse, csv, hashlib, json, pickle
from pathlib import Path
import numpy as np

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def endpoints(path):
    out=set()
    with open(path,encoding='utf-8',newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            out.update(r[k] for k in ('anchor_residue_id','true_target_residue_id','negative_residue_id'))
    return out
def canonical_array_sha(ids, emb):
    h=hashlib.sha256()
    for rid in ids:
        h.update(rid.encode('utf-8')); h.update(b'\0')
    h.update(np.asarray(emb,dtype='<f4',order='C').tobytes(order='C'))
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--model-id',required=True); ap.add_argument('--relations',required=True); ap.add_argument('--token-map',required=True); ap.add_argument('--residue-evidence'); ap.add_argument('--residue-dir'); ap.add_argument('--protoken-tokens'); ap.add_argument('--protoken-checkpoint'); ap.add_argument('--output',required=True); a=ap.parse_args()
    need=endpoints(a.relations); mapping={}
    witness={}
    if a.residue_evidence:
        with open(a.residue_evidence,encoding='utf-8',newline='') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                rid=r['residue_id']
                if rid in witness: raise ValueError('duplicate residue evidence: '+rid)
                witness[rid]=r['structure_object_id']
    with open(a.token_map,encoding='utf-8',newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['model_id']!=a.model_id or r['residue_id'] not in need: continue
            rid=r['residue_id']; obj=r['structure_object_id']; idx=int(r['token_index'])
            if rid in mapping: raise ValueError('duplicate model/residue mapping: '+rid)
            if idx < 0: raise ValueError('negative token index: '+rid)
            if rid.rpartition(':')[0] != obj: raise ValueError('residue/object identity mismatch: '+rid)
            if witness and witness.get(rid) != obj: raise ValueError('residue-evidence identity mismatch: '+rid)
            mapping[rid]=(obj,idx)
    if set(mapping)!=need: raise ValueError(f'planned endpoint mapping failure: {len(mapping)}/{len(need)}')
    ids=sorted(need); arrays={}
    if a.model_id=='ProToken':
        tok={str(k).removesuffix('_clean'):np.asarray(v,dtype=np.int64) for k,v in pickle.load(open(a.protoken_tokens,'rb')).items()}
        ck=pickle.load(open(a.protoken_checkpoint,'rb')); code=np.asarray(ck['params']['vq_tokenizer']['codebook'],dtype=np.float32)
        for rid in ids:
            o,i=mapping[rid]
            if o not in tok or i<0 or i>=len(tok[o]): raise ValueError('ProToken token closure failure: '+rid)
            arrays[rid]=code[tok[o][i]]
    else:
        root=Path(a.residue_dir); cache={}
        for rid in ids:
            o,i=mapping[rid]
            if o not in cache:
                with np.load(root/(o+'.npz'),allow_pickle=False) as z: cache[o]=np.asarray(z['residue_embedding'],dtype=np.float32)
            if i<0 or i>=len(cache[o]): raise ValueError('residue token closure failure: '+rid)
            arrays[rid]=cache[o][i]
    emb=np.asarray([arrays[x] for x in ids],dtype=np.float32)
    if emb.ndim!=2 or not np.isfinite(emb).all(): raise ValueError('endpoint embedding failure')
    Path(a.output).parent.mkdir(parents=True,exist_ok=True); np.savez_compressed(a.output,residue_ids=np.asarray(ids),embeddings=emb)
    rep={'status':'PASS','model_id':a.model_id,'endpoint_count':len(ids),'embedding_dimension':int(emb.shape[1]),'embeddings_sha256':sha(a.output),'canonical_array_sha256':canonical_array_sha(ids,emb),'mapping_qc':{'negative_indexes':0,'duplicate_model_residue_keys':0,'object_identity_mismatches':0,'residue_evidence_bound':bool(a.residue_evidence)}}
    Path(a.output).with_name('ENDPOINT_MATERIALIZATION_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n'); print(json.dumps(rep,sort_keys=True))
if __name__=='__main__': main()
