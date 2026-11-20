import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np, torch

def fasta(path):
    rows=[]; ident=None; seq=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if ident is not None: rows.append((ident,''.join(seq)))
            ident=line[1:].strip(); seq=[]
        else: seq.append(line.strip().upper())
    if ident is not None: rows.append((ident,''.join(seq)))
    return rows

p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--repo',required=True);p.add_argument('--out',required=True);p.add_argument('--device',default='cuda');a=p.parse_args()
sys.path.insert(0,a.repo)
from rinalmo.pretrained import get_pretrained_model
rows=fasta(a.fasta)
expected=[x[0] for x in rows]
t0=time.time()
model,alphabet=get_pretrained_model('giga-v1')
model.eval().to(a.device)
with torch.inference_mode(), torch.autocast(device_type="cuda", dtype=torch.float16):
    vals=[]
    for ident,seq in rows:
        tok=torch.tensor(alphabet.batch_tokenize([seq]),dtype=torch.long,device=a.device)
        h=model(tok)['representation']
        valid=h[0,1:len(seq)+1]
        if valid.shape[0] != len(seq): raise ValueError(f'coordinate mismatch {ident}')
        vals.append(valid.float().mean(0).cpu().numpy())
z=np.asarray(vals,dtype=np.float32)
if not np.isfinite(z).all(): raise ValueError('nonfinite')
np.savez(a.out,sequence_sha256=np.array(expected),vectors=z)
meta={'model_id':'RiNALMo-Giga','loader':'official_get_pretrained_model_giga-v1_strict','records':len(rows),'input_id_sha256':hashlib.sha256(('\n'.join(expected)+'\n').encode()).hexdigest(),'embedding_dim':int(z.shape[1]),'pooling':'mean valid nucleotide states; CLS/EOS excluded','dtype':'CUDA FP16 autocast; float32 output','device':a.device,'elapsed_seconds':time.time()-t0,'repo_commit':__import__('subprocess').check_output(['git','-C',a.repo,'rev-parse','HEAD'],text=True).strip(),'weight_sha256':hashlib.sha256((Path.home()/'.cache/rinalmo_pretrained/giga-v1.pt').read_bytes()).hexdigest()}
Path(str(a.out)+'.metrics.json').write_text(json.dumps(meta,sort_keys=True)+'\n')
