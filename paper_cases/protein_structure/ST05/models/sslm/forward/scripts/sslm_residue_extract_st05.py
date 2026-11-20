#!/usr/bin/env python3
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np, torch
from transformers import EsmConfig,EsmModel,EsmTokenizer
def fasta(p):
    d={}; k=None
    for x in p.read_text().splitlines():
        if x.startswith('>'): k=x[1:].split()[0]; d[k]=''
        elif x.strip(): d[k]+=x.strip()
    return d
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--fasta',type=Path,required=True); ap.add_argument('--model-dir',type=Path,required=True); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args()
    seq=fasta(a.fasta); a.output_dir.mkdir(parents=True,exist_ok=True); dev=torch.device('cuda:0'); tok=EsmTokenizer.from_pretrained(a.model_dir); model=EsmModel(EsmConfig.from_pretrained(a.model_dir))
    raw=torch.load(a.model_dir/'pytorch_model.bin',map_location='cpu',weights_only=False); pre='encoder.esm.'; state={k[len(pre):]:v for k,v in raw.items() if k.startswith(pre)}; loaded=model.load_state_dict(state,strict=False)
    bad=[k for k in loaded.missing_keys if not (k.startswith('pooler.') or k.startswith('contact_head.'))]; unexpected=[k for k in loaded.unexpected_keys if k!='embeddings.position_embeddings.weight']
    if bad or unexpected: raise RuntimeError('checkpoint closure failure')
    model=model.to(dev).eval(); records=[]; t=time.time()
    with torch.inference_mode(),torch.autocast('cuda',dtype=torch.float16):
        for i,(oid,s) in enumerate(sorted(seq.items()),1):
            dst=a.output_dir/(oid+'.npz')
            if dst.exists():
                with np.load(dst) as z: r=np.asarray(z['residue_embedding'])
            else:
                q=tok(s,return_tensors='pt',add_special_tokens=True).to(dev); h=model(**q).last_hidden_state[0]; mask=q['attention_mask'][0].bool(); mask[0]=False; mask[mask.nonzero()[-1]]=False; r=h[mask].float().cpu().numpy().astype(np.float32); np.savez_compressed(dst,residue_embedding=r)
            if r.shape[0]!=len(s) or not np.isfinite(r).all(): raise RuntimeError('3Di embedding closure failure')
            records.append({'object_id':oid,'residue_count':len(r),'npz_sha256':sha(dst)})
            if i%25==0: print(json.dumps({'done':i,'total':len(seq)}),flush=True)
    rep={'status':'PASS','model_interface_release_id':'SSLM_3DI_v1.0.0_20261010','objects':len(records),'pooling':'none; 3Di residue embeddings retained','checkpoint_sha256':sha(a.model_dir/'pytorch_model.bin'),'elapsed_seconds':time.time()-t,'records':records}; (a.output_dir/'FORWARD_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n')
if __name__=='__main__': main()
