#!/usr/bin/env python3
import argparse,csv,gzip,hashlib,json,math,os,re,sqlite3,subprocess,tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
COMMON={'matrix':'EBLOSUM62','gapopen':'10.0','gapextend':'0.5','endweight':'Y','endopen':'10.0','endextend':'0.5'}
AA=re.compile(r'^[ACDEFGHIKLMNPQRSTVWYBXZJUO]+$')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def readfa(p):
 d={};k=None;buf=[]
 for z in open(p):
  z=z.strip()
  if not z:continue
  if z.startswith('>'):
   if k is not None:
    if k in d:raise ValueError('duplicate FASTA id')
    d[k]=''.join(buf).upper()
   k=z[1:].split()[0];buf=[]
  else:buf.append(z)
 if k is not None:
  if k in d:raise ValueError('duplicate FASTA id')
  d[k]=''.join(buf).upper()
 if not d or any(not AA.fullmatch(s) for s in d.values()):raise ValueError('invalid protein FASTA')
 return d
def onefa(p,s):p.write_text('>x\n'+s+'\n')
def score(binary,a,b):
 with tempfile.TemporaryDirectory(prefix='ps05_nw_') as td:
  td=Path(td);fa,fb,o=td/'a.fa',td/'b.fa',td/'score.txt';onefa(fa,a);onefa(fb,b)
  cmd=[str(binary),'-asequence',str(fa),'-bsequence',str(fb),'-datafile',COMMON['matrix'],'-gapopen',COMMON['gapopen'],'-gapextend',COMMON['gapextend'],'-endweight',COMMON['endweight'],'-endopen',COMMON['endopen'],'-endextend',COMMON['endextend'],'-aformat','score','-outfile',str(o),'-auto']
  env=os.environ.copy();env['PATH']=str(binary.parent)+os.pathsep+env.get('PATH','');env.setdefault('EMBOSS_ACDROOT',str(binary.parent.parent/'share/EMBOSS/acd'));env.setdefault('EMBOSS_DATA',str(binary.parent.parent/'share/EMBOSS/data'))
  q=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
  if q.returncode:raise RuntimeError(q.stderr[:500])
  m=re.search(r'\(([-+0-9.eE]+)\)',o.read_text()) or re.search(r'Score:\s*([-+0-9.eE]+)',o.read_text())
  if not m:raise ValueError('unparseable EMBOSS score')
  return float(m.group(1))
def norm(raw,sa,sb):
 d=sa+sb
 if not math.isfinite(d) or d==0:raise ValueError('invalid self score')
 return 2*raw/d
def main():
 p=argparse.ArgumentParser();p.add_argument('--relations',type=Path,required=True);p.add_argument('--plan',type=Path,required=True);p.add_argument('--fasta',type=Path,required=True);p.add_argument('--tool',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--workers',type=int,default=32);a=p.parse_args()
 plan=json.load(open(a.plan));fields=plan['identity_fields'];ids=plan['leaf_identities']
 if len(ids)!=1628 or len({tuple(x) for x in ids})!=1628:raise ValueError('planned population cardinality drift')
 rows=[dict(zip(fields,x)) for x in ids];seq=readfa(a.fasta)
 rel={(r['anchor_domain_id'],r['positive_domain_id'],r['negative_domain_id']) for r in csv.DictReader(open(a.relations),delimiter='\t')}
 triples={(r['anchor_id'],r['positive_endpoint_id'],r['negative_endpoint_id']) for r in rows}
 if triples!=rel or len(rel)!=1628:raise ValueError('planned relation identity does not equal frozen relation ledger')
 endpoints=sorted({v for r in rows for v in (r['anchor_id'],r['positive_endpoint_id'],r['negative_endpoint_id'])})
 if len(endpoints)!=321 or set(endpoints)!=set(seq):raise ValueError('endpoint population drift')
 a.outdir.mkdir(parents=True,exist_ok=False); db=sqlite3.connect(a.outdir/'pair_index.sqlite');db.execute('CREATE TABLE pair(a TEXT,b TEXT,PRIMARY KEY(a,b))')
 for r in rows:
  for x,y in ((r['anchor_id'],r['positive_endpoint_id']),(r['anchor_id'],r['negative_endpoint_id'])):db.execute('INSERT OR IGNORE INTO pair VALUES (?,?)',tuple(sorted((x,y))))
 db.commit();pairs=list(db.execute('SELECT a,b FROM pair ORDER BY a,b'))
 with ThreadPoolExecutor(max_workers=a.workers) as ex:selfs=dict(ex.map(lambda x:(x,score(a.tool,seq[x],seq[x])),endpoints))
 with ThreadPoolExecutor(max_workers=a.workers) as ex:raw=dict(ex.map(lambda x:(x,score(a.tool,seq[x[0]],seq[x[1]])),pairs))
 with gzip.open(a.outdir/'PAIR_SCORE_CACHE.tsv.gz','wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['id_a','id_b','alignment_score','self_a','self_b','normalized_similarity'],delimiter='\t',lineterminator='\n');w.writeheader()
  for x,y in pairs:w.writerow({'id_a':x,'id_b':y,'alignment_score':repr(raw[(x,y)]),'self_a':repr(selfs[x]),'self_b':repr(selfs[y]),'normalized_similarity':repr(norm(raw[(x,y)],selfs[x],selfs[y]))})
 prov={'assay_id':plan['assay_id'],'planned_bundle_ref':'PS05_ASSAY_REFSET_RELEASE_v1.0.1','interface_coordinate_ref':'PS05_FORMAL_ENDPOINTS.faa','execution_attempt_ref':'PS05_NW_BLOSUM62_FORMAL_v1.0.0','measurement_contract_ref':plan['measurement_contract_ref']}
 leaf=a.outdir/'FORMAL_LEAF_SCORES.tsv';detail=a.outdir/'CONTRAST_SCORE_LEDGER.tsv.gz';wins=ties=losses=0
 with open(leaf,'w',newline='') as lf,gzip.open(detail,'wt',newline='') as df:
  fw=csv.DictWriter(lf,fieldnames=fields+list(prov)+['positive_score','negative_score'],delimiter='\t',lineterminator='\n');fw.writeheader();dw=csv.DictWriter(df,fieldnames=['anchor_id','positive_id','negative_id','S_anchor_positive','S_anchor_negative','delta','F_credit'],delimiter='\t',lineterminator='\n');dw.writeheader()
  for r in rows:
   ap=norm(raw[tuple(sorted((r['anchor_id'],r['positive_endpoint_id'])))],selfs[r['anchor_id']],selfs[r['positive_endpoint_id']]);an=norm(raw[tuple(sorted((r['anchor_id'],r['negative_endpoint_id'])))],selfs[r['anchor_id']],selfs[r['negative_endpoint_id']]);d=ap-an;f_credit=1. if d>0 else -1. if d<0 else 0.;wins+=f_credit==1;ties+=f_credit==0;losses+=f_credit==-1;fw.writerow({**r,**prov,'positive_score':repr(ap),'negative_score':repr(an)});dw.writerow({'anchor_id':r['anchor_id'],'positive_id':r['positive_endpoint_id'],'negative_id':r['negative_endpoint_id'],'S_anchor_positive':repr(ap),'S_anchor_negative':repr(an),'delta':repr(d),'F_credit':repr(f_credit)})
 summary={'kind':'PS05_NW_BLOSUM62_MEASUREMENT','status':'PASS','formal_result':False,'relation_count':len(rows),'endpoint_count':len(endpoints),'unique_pair_count':len(pairs),'wins':wins,'ties':ties,'losses':losses,'tool_sha256':sha(a.tool),'relation_ledger_sha256':sha(a.relations),'endpoint_fasta_sha256':sha(a.fasta),'planned_population_sha256':sha(a.plan),'pair_score_cache_sha256':sha(a.outdir/'PAIR_SCORE_CACHE.tsv.gz'),'leaf_ledger_sha256':sha(leaf),'detail_ledger_sha256':sha(detail),'parameters':COMMON,'fixed_prior':'EBLOSUM62','terminal_gaps_penalized':True}
 (a.outdir/'COMPARATOR_LEDGER_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':main()
