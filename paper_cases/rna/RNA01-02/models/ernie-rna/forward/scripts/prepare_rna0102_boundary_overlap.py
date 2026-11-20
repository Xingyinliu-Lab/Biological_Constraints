from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

def rows(p):
 r=[];k=None;s=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if k is not None:r.append((k,''.join(s)))
   k=x[1:].strip();s=[]
  else:s.append(x.strip().upper())
 if k is not None:r.append((k,''.join(s)))
 return r

def main():
 a=argparse.ArgumentParser();a.add_argument('--fasta',required=True);a.add_argument('--out',required=True);a.add_argument('--summary',required=True);a.add_argument('--window-nt',type=int,default=1022);z=a.parse_args()
 q=[x for x in rows(z.fasta) if len(x[1])>z.window_nt]
 with open(z.out,'w') as f:
  for i,s in q:f.write(f'>{i}\n{s}\n')
 Path(z.summary).write_text(json.dumps({'records':len(q),'window_nt':z.window_nt,'lengths':sorted(set(map(lambda x:len(x[1]),q))),'policy':'overlap merge uses windows [0,W) and [L-W,L) for every boundary sequence'},sort_keys=True)+'\n')
if __name__=='__main__':main()