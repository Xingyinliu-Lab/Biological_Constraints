#!/usr/bin/env python3
import argparse, csv, hashlib, json, os, shutil, sys, tempfile, time
from collections import defaultdict
from pathlib import Path
import numpy as np
import torch
from transformers import AutoModelForMaskedLM

WINDOW_BP = 1000

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def fasta(path):
    d={}; key=None
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            if key in d: raise ValueError('duplicate FASTA id')
            d[key]=''
        else: d[key]+=line.strip().upper()
    return d

def encode(model,ids,seqs):
    tasks=[]; weights=defaultdict(int); sums={}
    for oid in ids:
        s=seqs[oid]
        for i in range(0,len(s),WINDOW_BP):
            w=s[i:i+WINDOW_BP]; tasks.append((len(w),oid,w)); weights[oid]+=len(w)
    tasks.sort(key=lambda x:(x[0],x[1],x[2]))
    for start in range(0,len(tasks),8):
        batch=tasks[start:start+8]
        base_to_id={b:i for i,b in enumerate('ACGTN')}
        max_len=max(len(z[2]) for z in batch)
        encoded=[[base_to_id[b] for b in z[2]]+[4]*(max_len-len(z[2])) for z in batch]
        input_ids=torch.tensor(encoded,dtype=torch.long,device='cuda')
        with torch.inference_mode(),torch.autocast(device_type='cuda',dtype=torch.float16):
            hidden,_=model.backbone(input_ids)
        mask=torch.tensor([[True]*len(z[2])+[False]*(max_len-len(z[2])) for z in batch],device='cuda')
        v=(hidden.float()*mask.unsqueeze(-1)).sum(1)/mask.sum(1,keepdim=True).clamp_min(1)
        for (bp,oid,_),vec in zip(batch,v.cpu()): sums[oid]=sums.get(oid,0)+vec*bp
        if start==0 or start+len(batch)==len(tasks) or (start//8)%25==0:
            print(json.dumps({'model':'crossdna_8m','windows_done':min(start+len(batch),len(tasks)),'windows_total':len(tasks)}),flush=True)
    out=torch.stack([(sums[oid]/weights[oid]).float() for oid in ids]).numpy()
    if not np.isfinite(out).all(): raise ValueError('non-finite embeddings')
    return out,len(tasks)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--object-table',type=Path,required=True); ap.add_argument('--fasta',type=Path,required=True); ap.add_argument('--model-dir',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    if a.output.exists(): raise SystemExit('output path must not already exist')
    seqs=fasta(a.fasta)
    with a.object_table.open(encoding='utf-8',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    ids=[]
    for r in rows:
        oid=r['object_id']; s=seqs.get(oid)
        if s is None or len(s)!=int(r['length_bp']) or hashlib.sha256(s.encode()).hexdigest()!=r['sequence_sha256']: raise ValueError('object identity mismatch')
        ids.append(oid)
    model=AutoModelForMaskedLM.from_pretrained(a.model_dir,trust_remote_code=True,local_files_only=True).eval().half().cuda()
    for p in model.parameters(): p.requires_grad_(False)
    model.backbone.pretrain=False; model.backbone.for_representation=True
    started=time.time(); primary,nw=encode(model,ids,seqs); repeat,nr=encode(model,ids,seqs)
    a.output.parent.mkdir(parents=True,exist_ok=True); stage=Path(tempfile.mkdtemp(prefix=f'.{a.output.name}.staging.',dir=a.output.parent))
    try:
        np.savez_compressed(stage/'EMBEDDINGS.npz',object_ids=np.asarray(ids),vectors=primary)
        np.savez_compressed(stage/'REPEAT_EMBEDDINGS.npz',object_ids=np.asarray(ids),vectors=repeat)
        interface={'artifact':'DNA01_MODEL_INTERFACE','status':'EXECUTED_FROZEN_INTERFACE','model_key':'crossdna_8m','model_id':'CrossDNA_8.1M','representation':'backbone for_representation fused token matrix; token mean within each nonoverlapping 1000-bp window; bp-weighted window mean','window_bp':1000,'precision':'model float16; pooling and released vectors float32','operator_branch':'EMPIRICAL_SHARED_AFFINE','source_sha256':{'DNA01_OBJECT_TABLE.tsv':sha(a.object_table),'SOURCE_FASTA':sha(a.fasta)},'runtime':{'python':sys.version.split()[0],'torch':torch.__version__,'numpy':np.__version__}}
        (stage/'MODEL_INTERFACE.json').write_text(json.dumps(interface,indent=2,sort_keys=True)+'\n')
        validation={'status':'PASS','objects':len(ids),'dimension':int(primary.shape[1]),'windows':nw,'repeat_windows':nr,'repeat_max_abs_difference':float(np.max(np.abs(primary-repeat))),'all_finite':True,'elapsed_seconds':round(time.time()-started,3)}
        (stage/'FORWARD_VALIDATION.json').write_text(json.dumps(validation,indent=2,sort_keys=True)+'\n')
        names=['EMBEDDINGS.npz','REPEAT_EMBEDDINGS.npz','MODEL_INTERFACE.json','FORWARD_VALIDATION.json']; files={n:{'sha256':sha(stage/n),'bytes':(stage/n).stat().st_size} for n in names}
        (stage/'OUTPUT_MANIFEST.json').write_text(json.dumps({'status':'PASS','files':files},indent=2,sort_keys=True)+'\n')
        names.append('OUTPUT_MANIFEST.json'); (stage/'SHA256SUMS').write_text(''.join(f'{sha(stage/n)}  {n}\n' for n in names)); (stage/'RUN_COMPLETE').write_text('PASS\n')
        os.replace(stage,a.output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True); raise
    print(json.dumps(validation,sort_keys=True))

if __name__=='__main__': main()
