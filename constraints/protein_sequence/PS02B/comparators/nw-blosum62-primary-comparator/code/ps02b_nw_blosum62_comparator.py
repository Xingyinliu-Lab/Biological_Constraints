#!/usr/bin/env python3
"""Deterministic PS02B raw-sequence Needle+BLOSUM62 score-ledger producer."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math,os,re,subprocess,tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
P={'matrix':'EBLOSUM62','gapopen':'10.0','gapextend':'0.5','endweight':'Y','endopen':'10.0','endextend':'0.5'}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def op(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p,'rt')
def fasta(p):
 d={}; k=None; b=[]
 with op(p) as f:
  for z in f:
   z=z.strip()
   if not z: continue
   if z.startswith('>'):
    if k is not None:
     if k in d: raise ValueError('duplicate fasta id')
     d[k]=''.join(b).upper()
    k=z[1:].split()[0]; b=[]
   else:b.append(z)
 if k is not None:
  if k in d:raise ValueError('duplicate fasta id')
  d[k]=''.join(b).upper()
 if not d or any(not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+',v) for v in d.values()):raise ValueError('invalid protein FASTA')
 return d
def writefa(p,s):p.write_text('>x\n'+s+'\n')
def parse(t):
 m=re.search(r'\(([-+0-9.eE]+)\)',t) or re.search(r'Score:\s*([-+0-9.eE]+)',t)
 if not m:raise ValueError('unparseable needle score')
 return float(m.group(1))
def needle(binary,a,b):
 with tempfile.TemporaryDirectory(prefix='ps02b_nw_') as z:
  t=Path(z); fa,fb,o=t/'a.fa',t/'b.fa',t/'o.score';writefa(fa,a);writefa(fb,b)
  cmd=[str(binary),'-asequence',str(fa),'-bsequence',str(fb),'-datafile',P['matrix'],'-gapopen',P['gapopen'],'-gapextend',P['gapextend'],'-endweight',P['endweight'],'-endopen',P['endopen'],'-endextend',P['endextend'],'-aformat','score','-outfile',str(o),'-auto']
  e=os.environ.copy();e['PATH']=str(binary.parent)+os.pathsep+e.get('PATH','');e.setdefault('EMBOSS_ACDROOT',str(binary.parent.parent/'share/EMBOSS/acd'));e.setdefault('EMBOSS_DATA',str(binary.parent.parent/'share/EMBOSS/data'))
  q=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=e)
  if q.returncode:raise RuntimeError(q.stderr[:500])
  return parse(o.read_text())
def norm(s,aa,bb):
 d=aa+bb
 if not math.isfinite(d) or d==0:raise ValueError('invalid self score')
 return 2*s/d
def main():
 a=argparse.ArgumentParser();a.add_argument('--bundles',type=Path,required=True);a.add_argument('--fasta',type=Path,required=True);a.add_argument('--needle',type=Path,required=True);a.add_argument('--outdir',type=Path,required=True);a.add_argument('--workers',type=int,default=1);a.add_argument('--limit',type=int);x=a.parse_args();x.outdir.mkdir(parents=True,exist_ok=True)
 if x.workers<1:raise ValueError('workers must be positive')
 rows=list(csv.DictReader(op(x.bundles),delimiter='\t')); rows=rows if x.limit is None else rows[:x.limit]
 if not rows or not {'bundle_id','anchor','positive','negative'}<=set(rows[0]):raise ValueError('invalid bundle ledger')
 d=fasta(x.fasta); ids=sorted({u for r in rows for u in (r['anchor'],r['positive'],r['negative'])}); missing=set(ids)-set(d)
 if missing:raise ValueError(f'missing formal endpoint sequences: {len(missing)}')
 pairs=sorted({tuple(sorted((r['anchor'],q))) for r in rows for q in (r['positive'],r['negative'])})
 jobs=[(i,i) for i in ids]+pairs
 def job(k):return k,needle(x.needle,d[k[0]],d[k[1]])
 with ThreadPoolExecutor(max_workers=x.workers) as ex: scores=dict(ex.map(job,jobs))
 self={i:scores[(i,i)] for i in ids}; ps={z:scores[z] for z in pairs}
 cache=x.outdir/'PAIR_SCORE_CACHE.tsv.gz'
 with gzip.open(cache,'wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['id_a','id_b','needle_score','self_a','self_b','normalized_similarity'],delimiter='\t',lineterminator='\n');w.writeheader()
  for u,v in pairs:w.writerow({'id_a':u,'id_b':v,'needle_score':repr(ps[(u,v)]),'self_a':repr(self[u]),'self_b':repr(self[v]),'normalized_similarity':repr(norm(ps[(u,v)],self[u],self[v]))})
 leaf=x.outdir/'FORMAL_LEAF_SCORES.tsv.gz'; detail=x.outdir/'CONTRAST_SCORE_LEDGER.tsv.gz'; credits=[]
 with gzip.open(leaf,'wt',newline='') as f, gzip.open(detail,'wt',newline='') as g:
  lw=csv.DictWriter(f,fieldnames=['bundle_id','score_anchor_positive','score_anchor_negative'],delimiter='\t',lineterminator='\n');lw.writeheader();dw=csv.DictWriter(g,fieldnames=['bundle_id','S_anchor_positive','S_anchor_negative','delta','F_credit'],delimiter='\t',lineterminator='\n');dw.writeheader()
  for r in rows:
   ap=norm(ps[tuple(sorted((r['anchor'],r['positive'])))],self[r['anchor']],self[r['positive']]); an=norm(ps[tuple(sorted((r['anchor'],r['negative'])))],self[r['anchor']],self[r['negative']]); de=ap-an; cr=1. if de>0 else -1. if de<0 else 0.
   lw.writerow({'bundle_id':r['bundle_id'],'score_anchor_positive':repr(ap),'score_anchor_negative':repr(an)});dw.writerow({'bundle_id':r['bundle_id'],'S_anchor_positive':repr(ap),'S_anchor_negative':repr(an),'delta':repr(de),'F_credit':repr(cr)});credits.append(cr)
 out={'artifact_type':'PS02B_NW_BLOSUM62_COMPARATOR_SCORE_LEDGER','status':'PASS','formal_result':False,'bundle_count':len(rows),'endpoint_count':len(ids),'pair_count':len(ps),'worker_count':x.workers,'wins':credits.count(1.),'ties':credits.count(0.),'losses':credits.count(-1.),'needle_sha256':sha(x.needle),'fixed_method_prior':'EBLOSUM62 substitution matrix','params':P,'bundle_sha256':sha(x.bundles),'fasta_sha256':sha(x.fasta),'pair_score_cache_sha256':sha(cache),'leaf_ledger_sha256':sha(leaf),'detail_ledger_sha256':sha(detail)}
 (x.outdir/'COMPARATOR_LEDGER_SUMMARY.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
