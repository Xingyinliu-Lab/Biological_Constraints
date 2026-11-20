#!/usr/bin/env python3
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np, torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def coords(path):
    residues={}; order=[]
    for line in path.open(errors='replace'):
        if not line.startswith('ATOM'): continue
        atom=line[12:16].strip()
        if atom not in ('N','CA','C') or line[16] not in (' ','A'): continue
        key=(line[21],line[22:26],line[26])
        if key not in residues: residues[key]={}; order.append(key)
        residues[key].setdefault(atom,[float(line[30:38]),float(line[38:46]),float(line[46:54])])
    good=[np.asarray([residues[k][a] for a in ('N','CA','C')],dtype=np.float32) for k in order if all(a in residues[k] for a in ('N','CA','C'))]
    if not good: raise ValueError('no complete backbone residues')
    return np.stack(good)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',type=Path,required=True); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True); pdbs=sorted(a.input_dir.glob('*.pdb')); dev=torch.device('cuda:0'); torch.cuda.set_device(dev)
    model=ESM3.from_pretrained('esm3_sm_open_v1',device=dev).eval(); records=[]; t=time.time()
    with torch.inference_mode(), torch.autocast('cuda',dtype=torch.bfloat16):
        for i,p in enumerate(pdbs,1):
            dst=a.output_dir/(p.stem+'.npz'); xyz=coords(p)
            if dst.exists():
                with np.load(dst) as q: h=np.asarray(q['residue_embedding'])
            else:
                enc=model.encode(ESMProtein(sequence=None,coordinates=torch.from_numpy(xyz))); st=enc.structure.unsqueeze(0).to(dev); sc=enc.coordinates.unsqueeze(0).to(dev)
                seq=torch.full_like(st,model.tokenizers.sequence.mask_token_id); seq[:,0]=model.tokenizers.sequence.bos_token_id; seq[:,-1]=model.tokenizers.sequence.eos_token_id
                if not torch.all(seq[:,1:-1]==model.tokenizers.sequence.mask_token_id): raise RuntimeError('sequence masking failure')
                h=model.forward(sequence_tokens=seq,structure_tokens=st,structure_coords=sc).embeddings[0,1:-1].float().cpu().numpy().astype(np.float32)
                np.savez_compressed(dst,residue_embedding=h)
            if h.shape[0]!=len(xyz) or not np.isfinite(h).all(): raise RuntimeError('embedding closure failure')
            records.append({'object_id':p.stem,'residue_count':len(h),'npz_sha256':sha(dst)})
            if i%10==0: print(json.dumps({'done':i,'total':len(pdbs)}),flush=True)
    rep={'status':'PASS','model_interface_release_id':'ESM3_STRUCT_MODEL_INTERFACE_v1.0.0_FROZEN_20261011','objects':len(records),'native_sequence_tokens':0,'all_interior_sequence_tokens_masked':True,'pooling':'none; residue embeddings retained','elapsed_seconds':time.time()-t,'records':records}
    (a.output_dir/'FORWARD_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n')
if __name__=='__main__': main()
