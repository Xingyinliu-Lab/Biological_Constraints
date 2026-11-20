import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
import torch

def fasta(p):
    out=[]; ident=None; buf=[]
    for x in Path(p).read_text().splitlines():
        if x.startswith('>'):
            if ident is not None: out.append((ident,''.join(buf).upper()))
            ident=x[1:].strip(); buf=[]
        else: buf.append(x.strip())
    if ident is not None: out.append((ident,''.join(buf).upper()))
    if not out or len({x[0] for x in out})!=len(out): raise ValueError('invalid or duplicate FASTA identity')
    return out

def layout(n,w=1022,stride=766):
    if n<1: raise ValueError('empty sequence')
    if n<=w:return [0]
    v=list(range(0,n-w+1,stride)); end=n-w
    if v[-1]!=end:v.append(end)
    return v

a=argparse.ArgumentParser();a.add_argument('--fasta',required=True);a.add_argument('--repo',required=True);a.add_argument('--out',required=True);a.add_argument('--device',default='cuda');a.add_argument('--batch-size',type=int,default=1);a=a.parse_args()
sys.path.insert(0,a.repo)
from rinalmo.pretrained import get_pretrained_model
rows=fasta(a.fasta); ids=[x[0] for x in rows]; t=time.time()
model,alphabet=get_pretrained_model('giga-v1');model.eval().to(a.device)
vals=[]; windowed=0
with torch.inference_mode(), torch.autocast(device_type='cuda',dtype=torch.float16):
 for ident,seq in rows:
  starts=layout(len(seq)); windowed+=len(starts)>1; chunks=[seq[s:s+1022] for s in starts]; acc=None; counts=None
  for j in range(0,len(chunks),a.batch_size):
   cc=chunks[j:j+a.batch_size]; tok=torch.tensor(alphabet.batch_tokenize(cc),dtype=torch.long,device=a.device); h=model(tok)['representation']
   for st,k,s in zip(starts[j:j+a.batch_size],range(len(cc)),cc):
    v=h[k,1:len(s)+1]
    if v.shape[0]!=len(s):raise ValueError('coordinate mismatch')
    if acc is None:acc=torch.zeros((len(seq),v.shape[1]),device=a.device);counts=torch.zeros(len(seq),device=a.device)
    acc[st:st+len(s)]+=v;counts[st:st+len(s)]+=1
  if (counts==0).any():raise ValueError('uncovered nucleotide')
  vals.append((acc/counts[:,None]).mean(0).float().cpu().numpy())
z=np.asarray(vals,dtype=np.float32)
if not np.isfinite(z).all():raise ValueError('nonfinite')
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
np.savez(a.out,sequence_sha256=np.asarray(ids),vectors=z)
meta={'model_id':'RiNALMo-Giga','loader':'official_get_pretrained_model_giga-v1','records':len(ids),'input_id_sha256':hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest(),'embedding_dim':int(z.shape[1]),'pooling':'overlapping 1022-nt coordinate windows; states averaged per nucleotide before equal-nucleotide whole-object mean; CLS/EOS excluded','window_nt':1022,'stride_nt':766,'windowed_nodes':windowed,'dtype':'CUDA FP16 autocast; float32 output','repo_commit':__import__('subprocess').check_output(['git','-C',a.repo,'rev-parse','HEAD'],text=True).strip(),'weight_sha256':hashlib.sha256((Path.home()/'.cache/rinalmo_pretrained/giga-v1.pt').read_bytes()).hexdigest(),'elapsed_seconds':time.time()-t}
Path(str(a.out)+'.metrics.json').write_text(json.dumps(meta,sort_keys=True)+'\n')
