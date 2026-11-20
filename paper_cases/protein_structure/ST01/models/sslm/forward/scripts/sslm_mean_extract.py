#!/usr/bin/env python3
"""SSLM frozen sequential-encoder extractor over Foldseek 3Di strings."""
from __future__ import annotations
import argparse, csv, hashlib, json, time
from pathlib import Path
import numpy as np
import torch
from transformers import EsmConfig, EsmModel, EsmTokenizer

def read_fasta(p: Path):
    out={}; key=None; buf=[]
    for line in p.read_text().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if key is not None: out[key]=''.join(buf)
            key=line[1:].split()[0]; buf=[]
        else: buf.append(line)
    if key is not None: out[key]=''.join(buf)
    return out

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--fasta',type=Path,required=True); ap.add_argument('--objects',type=Path,required=True); ap.add_argument('--model-dir',type=Path,required=True); ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args(); seq=read_fasta(a.fasta)
    with a.objects.open(encoding='utf-8', newline='') as handle:
        object_map={Path(row['coordinate_locator']).stem:row['object_id'] for row in csv.DictReader(handle, delimiter='\t')}
    if len(seq)!=28: raise ValueError(f'expected 28 3Di records, got {len(seq)}')
    dev=torch.device('cuda:0'); tok=EsmTokenizer.from_pretrained(a.model_dir)
    config=EsmConfig.from_pretrained(a.model_dir); model=EsmModel(config)
    raw=torch.load(a.model_dir/'pytorch_model.bin',map_location='cpu',weights_only=False)
    prefix='encoder.esm.'; state={k[len(prefix):]:v for k,v in raw.items() if k.startswith(prefix)}
    loaded=model.load_state_dict(state,strict=False)
    bad_missing=[k for k in loaded.missing_keys if not (k.startswith('pooler.') or k.startswith('contact_head.'))]
    allowed_unexpected={'embeddings.position_embeddings.weight'}
    bad_unexpected=[k for k in loaded.unexpected_keys if k not in allowed_unexpected]
    if bad_missing or bad_unexpected:
        raise RuntimeError(f'checkpoint closure failure: missing={bad_missing}, unexpected={bad_unexpected}')
    model=model.to(dev).eval()
    ids=[]; vec=[]; t0=time.time()
    with torch.inference_mode(), torch.autocast('cuda',dtype=torch.float16):
        for oid,s in sorted(seq.items()):
            q=tok(s,return_tensors='pt',add_special_tokens=True).to(dev)
            h=model(**q).last_hidden_state[0]
            mask=q['attention_mask'][0].bool(); mask[0]=False; mask[mask.nonzero()[-1]]=False
            r=h[mask].float().cpu()
            if r.shape[0]!=len(s) or not torch.isfinite(r).all(): raise RuntimeError('3Di embedding closure failure')
            stem=oid.removesuffix('_clean')
            if stem not in object_map: raise ValueError(f'3Di record is not a formal ST01 object: {oid}')
            ids.append(object_map[stem]); vec.append(r.mean(0).numpy().astype(np.float32))
    a.output.parent.mkdir(parents=True,exist_ok=True); np.savez(a.output,object_ids=np.asarray(ids),embeddings=np.asarray(vec,dtype=np.float32))
    rep={'status':'PASS','model_interface_release_id':'SSLM_3DI_v1.0.0_20261010','input_count':len(ids),'embedding_dim':int(len(vec[0])),
         'input_modality':'Foldseek_3Di_only','readout':'pretrained_SSLM_sequential_encoder_last_hidden_state','pooling':'unweighted_3Di_token_mean_excluding_special_tokens',
         'checkpoint_sha256':sha(a.model_dir/'pytorch_model.bin'),'checkpoint_tensor_count':len(raw),'loaded_encoder_tensor_count':len(state),'checkpoint_load_closure':'PASS',
         'input_3di_fasta_sha256':sha(a.fasta),'output_sha256':sha(a.output),'elapsed_seconds':time.time()-t0}
    a.output.with_name('FORWARD_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n'); print(json.dumps(rep,sort_keys=True))
if __name__=='__main__': main()
