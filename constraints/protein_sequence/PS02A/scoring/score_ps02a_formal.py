#!/usr/bin/env python3
"""PS02A frozen scorer: leaf measurements to ID-level and aggregate results."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
import numpy as np
from collections import defaultdict
from pathlib import Path

ORDER=['negative_id','negative_family','negative_clan','positive_id','positive_family','anchor','anchor_parent','point_unit','anchor_family','anchor_clan']
FORMAL_DRAWS=20000
FORMAL_SEED=20261006
FORMAL_DELTA_NUM=0.0
def read_tsv(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def mean(xs):return math.fsum(xs)/len(xs)
def aggregate(rows, ledger=None):
 cur=rows
 for i in range(len(ORDER)):
  keys=ORDER[i:]; groups=defaultdict(list)
  for r in cur:groups[tuple(r[k] for k in keys)].append(r['value'])
  cur=[]
  for k,v in groups.items():
   item={**dict(zip(keys,k)),'value':mean(v)}
   cur.append(item)
   if ledger is not None:ledger.append({'level':ORDER[i],'group_keys':dict(zip(keys,k)),'child_count':len(v),'weights':[1.0/len(v)]*len(v),'value':item['value']})
 result=mean([r['value'] for r in cur])
 if ledger is not None:ledger.append({'level':'OVERALL','group_keys':{},'child_count':len(cur),'weights':[1.0/len(cur)]*len(cur),'value':result})
 return result
def pct(xs,q):
 y=sorted(xs);pos=(len(y)-1)*q;lo=int(pos);hi=min(lo+1,len(y)-1);return y[lo]+(y[hi]-y[lo])*(pos-lo)
def write_individual_scores(path, rows):
 Path(path).parent.mkdir(parents=True,exist_ok=True)
 fields=['bundle_id',*ORDER,'score_anchor_positive','score_anchor_negative','score_difference','F_credit','bootstrap_unit']
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rows:w.writerow({k:r[k] for k in fields})
def main(a):
 bundles=read_tsv(a.bundles);points={r['bundle_id']:r['group_id'] for r in read_tsv(a.point_membership)};units={r['bundle_id']:r['anchor_clan_bootstrap_unit'] for r in read_tsv(a.uncertainty_membership)}; scores={r['bundle_id']:r for r in read_tsv(a.leaf_scores)}
 ids={r['bundle_id'] for r in bundles}
 if set(scores)!=ids or set(points)!=ids or set(units)!=ids:raise SystemExit('planned bundle membership mismatch')
 if any(units[b['bundle_id']]!=b['anchor_clan'] for b in bundles):raise SystemExit('uncertainty membership is not the bound PS02A v1.0.1 anchor-clan mapping')
 rows=[]
 for b in bundles:
  s=scores[b['bundle_id']];p=float(s['score_anchor_positive']);n=float(s['score_anchor_negative']);delta=FORMAL_DELTA_NUM
  if not(math.isfinite(p) and math.isfinite(n)):raise SystemExit('nonfinite measurement')
  value=1.0 if p-n>delta else -1.0 if p-n<-delta else 0.0
  rows.append({'negative_id':b['negative'],'negative_family':b['negative_family'],'negative_clan':b['negative_clan'],'positive_id':b['positive'],'positive_family':b['positive_family'],'anchor':b['anchor'],'anchor_parent':b['anchor_parent'],'point_unit':points[b['bundle_id']],'anchor_family':b['anchor_family'],'anchor_clan':b['anchor_clan'],'bootstrap_unit':units[b['bundle_id']],'value':value,'F_credit':value,'score_anchor_positive':p,'score_anchor_negative':n,'score_difference':p-n,'bundle_id':b['bundle_id']})
 clans=sorted({r['bootstrap_unit'] for r in rows});by_clan=defaultdict(list)
 for r in rows:by_clan[r['bootstrap_unit']].append(r)
 aggregation=[];point=aggregate(rows,aggregation);clan_score={c:aggregate(v) for c,v in by_clan.items()}
 if not aggregation or abs(aggregation[-1]['value']-point)>1e-12:raise SystemExit('dependency aggregation reconstruction drift')
 rng=np.random.Generator(np.random.PCG64(FORMAL_SEED));boot=[mean([clan_score[str(c)] for c in rng.choice(clans,size=len(clans),replace=True)]) for _ in range(FORMAL_DRAWS)]
 reuse=defaultdict(set)
 for r in rows:reuse[r['negative_id']].add(r['anchor_clan'])
 reused={n:c for n,c in reuse.items() if len(c)>1}; sens=[]
 for n,cs in reused.items():
  keep=[r for r in rows if r['negative_id']!=n];keep_clans={r['anchor_clan'] for r in keep}
  score=None if not keep else aggregate(keep)
  sens.append({'negative_id':n,'clans_using_negative':len(cs),'bundles_removed':len(rows)-len(keep),'retained_clans':len(keep_clans),'status':'NOT_ESTIMABLE_EMPTY_RETAINED_POPULATION' if score is None else 'PASS','F_without_negative':score,'delta_F':None if score is None else score-point})
 ci_f=[pct(boot,.025),pct(boot,.975)];valid=[x for x in sens if x['delta_F'] is not None];reuse_counts=[len(x) for x in reused.values()]
 Path(a.aggregation_ledger).parent.mkdir(parents=True,exist_ok=True)
 with gzip.open(a.aggregation_ledger,'wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['level','group_keys','child_count','weights','value'],delimiter='\t',lineterminator='\n');w.writeheader()
  for row in aggregation:w.writerow({**row,'group_keys':json.dumps(row['group_keys'],sort_keys=True,separators=(',',':')),'weights':json.dumps(row['weights'],separators=(',',':'))})
 aggregation_sha256=hashlib.sha256(Path(a.aggregation_ledger).read_bytes()).hexdigest()
 out={'formal_result':True,'formal_interface':'F_ONLY','score_transform_version':'PS02A_v1.0.2_F_ONLY','F':point,'n_bundles':len(rows),'n_anchor_clans':len(clans),'n_bootstrap_units':len(clans),'coverage':{'planned_bundles':len(bundles),'scored_bundles':len(rows),'fraction':1.0},'delta_num':FORMAL_DELTA_NUM,'bootstrap':{'rng':'PCG64','draws':FORMAL_DRAWS,'seed':FORMAL_SEED,'F_95_CI':ci_f},'dependency_aggregation_ledger':{'path':str(a.aggregation_ledger),'sha256':aggregation_sha256,'terminal_F':aggregation[-1]['value']},'cross_clan_negative_reuse':{'negative_count':len(reused),'reuse_clan_count_median':pct(reuse_counts,.5) if reuse_counts else None,'reuse_clan_count_max':max(reuse_counts) if reuse_counts else None,'sensitivity_delta_F_quantiles':[pct([x['delta_F'] for x in valid],q) for q in (.025,.5,.975)] if valid else None,'sensitivity_max_abs_delta_F':max((abs(x['delta_F']) for x in valid),default=None),'sensitivity':sens}}
 Path(a.output).parent.mkdir(parents=True,exist_ok=True);Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
 write_individual_scores(a.individual_scores,rows)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--bundles',required=True);p.add_argument('--point-membership',required=True);p.add_argument('--leaf-scores',required=True);p.add_argument('--uncertainty-membership',required=True);p.add_argument('--output',required=True);p.add_argument('--individual-scores',required=True);p.add_argument('--aggregation-ledger',required=True);main(p.parse_args())




