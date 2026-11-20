#!/usr/bin/env python3
import argparse,csv,hashlib,json,pickle
from pathlib import Path
import numpy as np

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--relations',type=Path,required=True); ap.add_argument('--tokens',type=Path,required=True); ap.add_argument('--checkpoint',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
 raw=pickle.load(a.tokens.open('rb')); tokens={str(k).removesuffix('_clean'):np.asarray(v,dtype=np.int64) for k,v in raw.items()}
 ck=pickle.load(a.checkpoint.open('rb')); code=np.asarray(ck['params']['vq_tokenizer']['codebook'],dtype=np.float32)
 with a.relations.open(encoding='utf-8',newline='') as h: rows=list(csv.DictReader(h,delimiter='\t'))
 bases=sorted({r['base_object_id'] for r in rows}); relids=[r['relation_id'] for r in rows]
 needed={x+'__BASE' for x in bases}|{r['variant_object_id'] for r in rows}
 if set(tokens)!=needed: raise ValueError(f'ProToken input closure failure missing={len(needed-set(tokens))} extra={len(set(tokens)-needed)}')
 def emb(k):
  idx=tokens[k]
  if idx.ndim!=1 or len(idx)==0 or np.any(idx<0) or np.any(idx>=len(code)): raise ValueError(f'invalid tokens {k}')
  return code[idx].mean(0,dtype=np.float64).astype(np.float32)
 B=np.asarray([emb(x+'__BASE') for x in bases],dtype=np.float32); V=np.asarray([emb(r['variant_object_id']) for r in rows],dtype=np.float32)
 if not np.isfinite(B).all() or not np.isfinite(V).all(): raise ValueError('nonfinite')
 a.output.parent.mkdir(parents=True,exist_ok=True); np.savez_compressed(a.output,base_object_ids=np.asarray(bases),base_embeddings=B,relation_ids=np.asarray(relids),variant_embeddings=V)
 rep={'status':'PASS','model_id':'ProToken','base_objects':len(bases),'relations':len(rows),'dimension':int(B.shape[1]),'checkpoint_sha256':sha(a.checkpoint),'token_artifact_sha256':sha(a.tokens),'embedding_sha256':sha(a.output),'pooling':'float64 mean of frozen codebook rows then float32'}; a.output.with_name('FORWARD_REPORT.json').write_text(json.dumps(rep,indent=2)+'\n'); print(json.dumps(rep,sort_keys=True))
if __name__=='__main__': main()
