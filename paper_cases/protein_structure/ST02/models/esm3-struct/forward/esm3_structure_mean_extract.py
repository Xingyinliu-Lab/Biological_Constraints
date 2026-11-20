#!/usr/bin/env python3
"""ESM3 structure-only extractor for the frozen ST02 structure population."""
from __future__ import annotations
import argparse, csv, hashlib, json, time
from pathlib import Path
import numpy as np
import torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein

INTERFACE_ID='ESM3_STRUCT_v1.0.0_20261010'

def sha(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def coords_from_pdb(path: Path) -> np.ndarray:
    residues={}; order=[]
    for line in path.open(errors='replace'):
        if not line.startswith('ATOM'): continue
        atom=line[12:16].strip()
        if atom not in ('N','CA','C') or line[16] not in (' ','A'): continue
        key=(line[21],line[22:26],line[26])
        if key not in residues: residues[key]={}; order.append(key)
        residues[key].setdefault(atom,[float(line[30:38]),float(line[38:46]),float(line[46:54])])
    good=[np.asarray([residues[k][a] for a in ('N','CA','C')],dtype=np.float32)
          for k in order if all(a in residues[k] for a in ('N','CA','C'))]
    if not good: raise ValueError('no complete backbone residues')
    return np.stack(good)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',type=Path,required=True); ap.add_argument('--objects',type=Path,required=True); ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args(); pdbs=sorted(a.input_dir.glob('*.pdb'))
    with a.objects.open(encoding='utf-8-sig', newline='') as handle: object_map={Path(row['coordinate_locator']).stem:row['record_id'] for row in csv.DictReader(handle,delimiter='\t')}
    if len(pdbs)!=509: raise ValueError(f'expected 509 inputs, got {len(pdbs)}')
    dev=torch.device('cuda:0'); torch.cuda.set_device(dev); torch.cuda.reset_peak_memory_stats(dev)
    model=ESM3.from_pretrained('esm3_sm_open_v1',device=dev).eval()
    ids=[]; vec=[]; audit=[]; t0=time.time()
    with torch.inference_mode(), torch.autocast('cuda',dtype=torch.bfloat16):
        for p in pdbs:
            xyz=torch.from_numpy(coords_from_pdb(p))
            enc=model.encode(ESMProtein(sequence=None,coordinates=xyz))
            st=enc.structure.unsqueeze(0).to(dev); sc=enc.coordinates.unsqueeze(0).to(dev)
            seq=torch.full_like(st,model.tokenizers.sequence.mask_token_id)
            seq[:,0]=model.tokenizers.sequence.bos_token_id; seq[:,-1]=model.tokenizers.sequence.eos_token_id
            if not torch.all(seq[:,1:-1] == model.tokenizers.sequence.mask_token_id): raise RuntimeError('sequence masking failure')
            out=model.forward(sequence_tokens=seq,structure_tokens=st,structure_coords=sc)
            h=out.embeddings[0,1:-1].float().cpu()
            if h.shape[0]!=xyz.shape[0] or not torch.isfinite(h).all(): raise RuntimeError('embedding closure failure')
            v=h.mean(0).numpy().astype(np.float32)
            if p.stem not in object_map: raise ValueError(f'coordinate file is not a formal ST02 object: {p.name}')
            ids.append(object_map[p.stem]); vec.append(v); audit.append({'residues':int(h.shape[0]),'interior_sequence_mask_fraction':1.0})
    a.output.parent.mkdir(parents=True,exist_ok=True)
    np.savez(a.output,object_ids=np.asarray(ids),embeddings=np.asarray(vec,dtype=np.float32))
    report={'status':'PASS','model_interface_release_id':INTERFACE_ID,'input_count':len(ids),'embedding_dim':int(len(vec[0])),
            'track_policy':'experimental_structure_only; native sequence absent; sequence interior all mask; function/SS/SASA omitted',
            'structure_representation':'ESM3 structure tokens plus backbone coordinates','pooling':'unweighted residue mean excluding boundaries',
            'compute_dtype':'bfloat16','output_dtype':'float32','input_token_audit':{'objects':len(audit),'all_interior_sequence_tokens_masked':True,'native_sequence_tokens':0},
            'output_sha256':sha(a.output),'elapsed_seconds':time.time()-t0,'peak_gpu_mib':torch.cuda.max_memory_allocated(dev)/2**20}
    a.output.with_name('FORWARD_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='input_token_audit'},sort_keys=True))
if __name__=='__main__': main()
