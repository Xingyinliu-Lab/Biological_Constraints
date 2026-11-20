#!/usr/bin/env python3
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import torch
from transformers import AutoModel

def fasta(path):
    ident=None; parts=[]
    for line in Path(path).read_text().splitlines():
        if not line: continue
        if line.startswith('>'):
            if ident is not None: yield ident,''.join(parts)
            ident=line[1:].split()[0]; parts=[]
        else: parts.append(line.strip().upper())
    if ident is not None: yield ident,''.join(parts)

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--fasta',required=True); a.add_argument('--output',required=True)
    a.add_argument('--mode',choices=('global','positions'),required=True)
    a.add_argument('--model-id',default='antichronology/orthrus-4-track')
    a.add_argument('--device',default='cuda'); a.add_argument('--batch-size',type=int,default=1)
    x=a.parse_args(); Path(x.output).parent.mkdir(parents=True, exist_ok=True); device=torch.device(x.device)
    records=list(fasta(x.fasta)); ids=[]
    for i,s in records:
        if not s or set(s)-set('ACGU'): raise ValueError(f'{i}: noncanonical RNA input')
        if hashlib.sha256(s.encode()).hexdigest()!=i: raise ValueError(f'{i}: sequence_sha256 mismatch')
        ids.append(i)
    if len(set(ids)) != len(ids): raise ValueError('duplicate sequence_sha256')
    model=AutoModel.from_pretrained(x.model_id,trust_remote_code=True,local_files_only=True).to(device).eval()
    if x.mode=='global':
        out=[]
        for start in range(0,len(records),x.batch_size):
            b=records[start:start+x.batch_size]
            oh=[model.seq_to_oh(s.replace('U','T')) for _,s in b]
            lengths=torch.tensor([z.shape[0] for z in oh],device=device)
            L=max(int(z.shape[0]) for z in oh)
            inp=torch.zeros((len(oh),L,4),dtype=oh[0].dtype,device=device)
            for k,z in enumerate(oh): inp[k,:z.shape[0]]=z.to(device)
            with torch.no_grad(): z=model.representation(inp,lengths,channel_last=True)
            out.append(z.detach().float().cpu().numpy())
        np.savez_compressed(x.output,sequence_sha256=np.array(ids),vectors=np.concatenate(out,axis=0))
    else:
        out={}
        for i,s in records:
            oh=model.seq_to_oh(s.replace('U','T')).unsqueeze(0).to(device)
            with torch.no_grad(): z=model.representation_unpooled(oh,channel_last=True)
            z=z.squeeze(0).detach().float().cpu().numpy()
            if z.shape[0]!=len(s): raise ValueError(f'{i}: nucleotide coordinate mismatch {z.shape[0]} != {len(s)}')
            out[i]=z
        np.savez_compressed(x.output,**out)
    meta={'model_id':x.model_id,'mode':x.mode,'sequence_count':len(records),'sequence_sha256_verified':True,'adapter':'RNA_U_TO_DNA_T','output':str(x.output)}
    Path(str(x.output)+'.metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__': main()
