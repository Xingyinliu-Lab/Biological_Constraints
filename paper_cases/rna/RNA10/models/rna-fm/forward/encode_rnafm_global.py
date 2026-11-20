from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
import torch

def fasta(path):
    rows=[]; ident=None; seq=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if ident is not None: rows.append((ident,''.join(seq)))
            ident=line[1:].strip(); seq=[]
        else: seq.append(line.strip().upper())
    if ident is not None: rows.append((ident,''.join(seq)))
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--weights',required=True);p.add_argument('--source-root',required=True);p.add_argument('--output',required=True);p.add_argument('--batch-size',type=int,default=4);p.add_argument('--device',default='cuda');a=p.parse_args()
    sys.path.insert(0,a.source_root);import fm
    rows=fasta(a.fasta)
    if not rows: raise ValueError('empty input')
    # RNA-FM reserves two positions for BOS/EOS.  For the frozen 1024-nt
    # population, chunks are only required for 9 boundary-length nodes.
    max_nt=1022
    model,alpha=fm.pretrained.rna_fm_t12(a.weights);model.eval().to(a.device);convert=alpha.get_batch_converter();vals=[];ids=[];t=time.time()
    def encode_batch(batch):
        _,_,tokens=convert(batch);h=model(tokens.to(a.device),repr_layers=[12])['representations'][12]
        return [h[j,1:len(seq)+1].float().cpu().numpy() for j,(_,seq) in enumerate(batch)]
    with torch.inference_mode():
      short=[]
      for ident,seq in rows:
        if len(seq)<=max_nt: short.append((ident,seq));continue
        parts=[seq[i:i+max_nt] for i in range(0,len(seq),max_nt)]
        states=[]
        for q in parts: states.extend(encode_batch([(ident,q)]))
        vals.append(np.concatenate(states,axis=0).mean(0,dtype=np.float32));ids.append(ident)
      for i in range(0,len(short),a.batch_size):
        b=short[i:i+a.batch_size]
        for (ident,_),h in zip(b,encode_batch(b)): vals.append(h.mean(0,dtype=np.float32));ids.append(ident)
        if i%(a.batch_size*100)==0 or i+len(b)==len(short):print(json.dumps({'short_done':i+len(b),'short_total':len(short),'windowed_boundary_nodes':len(rows)-len(short)}),flush=True)
    z=np.stack(vals).astype(np.float32)
    if not np.isfinite(z).all():raise ValueError('nonfinite')
    order=np.argsort(np.array(ids));z=z[order];ids=np.array(ids)[order]
    Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    np.savez(a.output,sequence_sha256=ids,vectors=z)
    m={'model_id':'RNA-FM','records':len(ids),'embedding_dim':int(z.shape[1]),'pooling':'mean valid nucleotide states; BOS/EOS excluded; 1024-nt boundary nodes tiled at 1022 nt before whole-object mean','windowed_boundary_nodes':len(rows)-len(short),'elapsed_seconds':time.time()-t,'weights_sha256':hashlib.sha256(Path(a.weights).read_bytes()).hexdigest(),'fasta_sha256':hashlib.sha256(Path(a.fasta).read_bytes()).hexdigest()}
    Path(str(a.output)+'.metrics.json').write_text(json.dumps(m,sort_keys=True)+'\n')
if __name__=='__main__':main()
