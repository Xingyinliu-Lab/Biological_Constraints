import hashlib,json
from pathlib import Path
import numpy as np
K=range(1,7); B=17; alphabet='ACGT'
def fasta(p):
 o=[]; i=None; b=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if i:o.append((i,''.join(b).upper().replace('U','T')))
   i=x[1:].strip();b=[]
  else:b.append(x.strip())
 if i:o.append((i,''.join(b).upper().replace('U','T')))
 return o
code={};off=0
for k in K:
 for b in range(B):
  for n in range(4**k):code[(k,b,n)]=off;off+=1
def enc(s):
 v=np.zeros(off,np.float32); L=len(s)
 for k in K:
  for p in range(L-k+1):
   q=s[p:p+k]
   if any(x not in alphabet for x in q):continue
   n=0
   for x in q:n=4*n+alphabet.index(x)
   b=min(B-1,(p*B)//max(1,L-k+1));v[code[(k,b,n)]]+=1
 return v
rows=fasta(Path(__import__('sys').argv[1])); ids=[x[0] for x in rows]
X=np.stack([enc(s) for _,s in rows]);output=Path(__import__('sys').argv[2]);output.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(output,sequence_sha256=np.array(ids),vectors=X)
Path(str(output)+'.metadata.json').write_text(json.dumps({'comparator_id':'CMP_RNA10_RAW_POSITIONAL_KMER_RIDGE_v0.1.0','input':'raw reporter sequence only; U-to-T lexical adapter','feature':'17 equal relative-position bins x k-mer counts for k=1..6','dimension':int(X.shape[1]),'records':len(ids),'feature_sha256':hashlib.sha256(X.tobytes()).hexdigest()},sort_keys=True)+'\n')
