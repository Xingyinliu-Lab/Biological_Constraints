import argparse,sys,json,time,hashlib
from pathlib import Path
import torch

def fa(p):
 o=[];i=None;b=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if i:o.append((i,''.join(b).upper()))
   i=x[1:].strip();b=[]
  else:b.append(x.strip())
 if i:o.append((i,''.join(b).upper()))
 return o
def lay(n):
 if n<=1022:return [0]
 v=list(range(0,n-1022+1,766));e=n-1022
 if v[-1]!=e:v.append(e)
 return v
a=argparse.ArgumentParser();a.add_argument('--fasta',required=True);a.add_argument('--repo',required=True);a.add_argument('--out',required=True);a.add_argument('--device',default='cuda');a=a.parse_args();sys.path.insert(0,a.repo)
from rinalmo.pretrained import get_pretrained_model
m,ab=get_pretrained_model('giga-v1');m.eval().to(a.device);out={};t=time.time()
with torch.inference_mode(),torch.autocast(device_type='cuda',dtype=torch.float16):
 for ident,s in fa(a.fasta):
  acc=cnt=None
  for st in lay(len(s)):
   c=s[st:st+1022]; h=m(torch.tensor(ab.batch_tokenize([c]),dtype=torch.long,device=a.device))['representation'][0,1:len(c)+1]
   if h.shape[0]!=len(c):raise ValueError('coordinate mismatch')
   if acc is None:acc=torch.zeros((len(s),h.shape[1]),device=a.device);cnt=torch.zeros(len(s),device=a.device)
   acc[st:st+len(c)]+=h;cnt[st:st+len(c)]+=1
  if (cnt==0).any():raise ValueError('uncovered')
  out[ident]=(acc/cnt[:,None]).float().cpu()
if not all(torch.isfinite(x).all() for x in out.values()):raise ValueError('nonfinite')
torch.save(out,a.out)
Path(str(a.out)+'.metrics.json').write_text(json.dumps({'model_id':'RiNALMo-Giga','loader':'official_get_pretrained_model_giga-v1','records':len(out),'representation':'per-nucleotide hidden state vectors','pooling':'1022-nt windows stride 766; coordinate overlap mean; CLS/EOS excluded','elapsed_seconds':time.time()-t},sort_keys=True)+'\n')
