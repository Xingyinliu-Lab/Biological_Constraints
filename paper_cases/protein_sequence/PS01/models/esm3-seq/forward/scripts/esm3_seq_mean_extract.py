#!/usr/bin/env python3
"""ESM3 sequence-only frozen-interface extractor; emits no biological text."""
from __future__ import annotations
import argparse, csv, hashlib, json, os, time
from pathlib import Path
import numpy as np
import torch
from safetensors.torch import save_file
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein

WINDOW=1023
STRIDE=767
INTERFACE_ID='ESM3_SEQ_v1.0.0_20261010'

def sha(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def read_fasta(p: Path):
    opener=__import__('gzip').open if p.suffix=='.gz' else open
    out={}; key=None; buf=[]
    with opener(p,'rt',encoding='utf-8') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if key is not None:
                    if key in out: raise RuntimeError('duplicate identifier')
                    out[key]=''.join(buf)
                key=line[1:].split()[0]; buf=[]
            else: buf.append(line.upper())
    if key is not None:
        if key in out: raise RuntimeError('duplicate identifier')
        out[key]=''.join(buf)
    if not out: raise RuntimeError('empty FASTA')
    return out

def spans(n: int):
    if n<=WINDOW: return [(0,n)]
    starts=list(range(0,n-WINDOW+1,STRIDE))
    if starts[-1] != n-WINDOW: starts.append(n-WINDOW)
    return [(s,min(n,s+WINDOW)) for s in starts]

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--fasta',required=True,type=Path)
    a.add_argument('--output-dir',required=True,type=Path)
    a.add_argument('--id-column',required=True)
    args=a.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    seq=read_fasta(args.fasta)
    dev=torch.device('cuda:0')
    if not torch.cuda.is_available(): raise RuntimeError('CUDA required')
    torch.cuda.set_device(dev); torch.cuda.reset_peak_memory_stats(dev)
    model=ESM3.from_pretrained('esm3_sm_open_v1',device=dev).eval()
    vectors={}; index=[]; tiles=0; t0=time.time()
    with torch.inference_mode(), torch.autocast('cuda',dtype=torch.bfloat16):
        for k,s in seq.items():
            acc=None; count=None
            for left,right in spans(len(s)):
                tok=model.encode(ESMProtein(sequence=s[left:right])).sequence.unsqueeze(0).to(dev)
                h=model.forward(sequence_tokens=tok).embeddings[0,1:-1].float().cpu()
                if h.shape[0] != right-left: raise RuntimeError('residue length mismatch')
                if acc is None:
                    acc=torch.zeros((len(s),h.shape[1]),dtype=torch.float32); count=torch.zeros(len(s),dtype=torch.float32)
                acc[left:right]+=h; count[left:right]+=1; tiles+=1
            if (count==0).any(): raise RuntimeError('uncovered residue')
            v=(acc/count[:,None]).mean(0).contiguous()
            if not torch.isfinite(v).all(): raise RuntimeError('non-finite embedding')
            ekey=f'embedding_{len(vectors)+1:06d}'; vectors[ekey]=v
            index.append((k,ekey,len(s),v.numel(),hashlib.sha256(v.numpy().tobytes()).hexdigest()))
    dest=args.output_dir/'embeddings.safetensors'; save_file(vectors,str(dest),metadata={
      'model_interface_release_id':INTERFACE_ID,'pooling':'residue_mean','sequence_window_residues':str(WINDOW),'sequence_stride_residues':str(STRIDE),'compute_dtype':'bfloat16','output_dtype':'float32'})
    ip=args.output_dir/'embedding_index.tsv'
    with ip.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow((args.id_column,'embedding_key','sequence_length','embedding_dim','embedding_sha256')); w.writerows(index)
    report={'status':'PASS','model_interface_release_id':INTERFACE_ID,'track_policy':'sequence_only; all non-sequence tracks omitted and model-default masked','pooling':'residue_mean','sequence_window_residues':WINDOW,'sequence_stride_residues':STRIDE,'compute_dtype':'bfloat16','output_dtype':'float32','input_count':len(index),'tile_count':tiles,'embedding_dim':index[0][3],'input_fasta_sha256':sha(args.fasta),'embedding_sha256':sha(dest),'embedding_index_sha256':sha(ip),'elapsed_seconds':time.time()-t0,'peak_gpu_mib':torch.cuda.max_memory_allocated(dev)/2**20}
    (args.output_dir/'FORWARD_REPORT.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
if __name__=='__main__': main()
