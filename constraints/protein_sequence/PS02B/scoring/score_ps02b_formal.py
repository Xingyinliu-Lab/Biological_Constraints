#!/usr/bin/env python3
"""Frozen PS02B hierarchy scorer for the L AssayRefSetRelease."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np

FORMAL_DRAWS=20000
FORMAL_SEED=20261006
FORMAL_DELTA_NUM=0.0

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def op(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p,'rt')
def rd(p):
 with op(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def mean(rows,keys):
 g=defaultdict(list)
 for r in rows:g[tuple(r[k] for k in keys)].append(r)
 out=[]
 for k,v in g.items():
  z=dict(v[0]);z.update(zip(keys,k));z['value']=math.fsum(float(x['value']) for x in v)/len(v);out.append(z)
 return out
def anchor_values(rows):
 # The sealed L release defines each role-local source unit as its source parent.
 for r in rows:
  r['negative_source_unit']=r['negative_parent'];r['positive_source_unit']=r['positive_parent'];r['anchor_source_unit']=r['anchor_parent']
 x=mean(rows,['anchor','positive','negative_hsf','negative'])
 x=mean(x,['anchor','positive','negative_hsf','negative_parent'])
 x=mean(x,['anchor','positive','negative_hsf','negative_source_unit'])
 x=mean(x,['anchor','positive','negative_hsf'])
 x=mean(x,['anchor','positive'])
 x=mean(x,['anchor','positive_parent'])
 x=mean(x,['anchor','positive_source_unit'])
 return {r['anchor']:r for r in mean(x,['anchor'])}
def hierarchy(meta,groups):
 cur=meta; comp=np.array(groups,dtype=np.int32); levels=[]
 for fields in [('topology','anchor_hsf','anchor_parent'),('topology','anchor_hsf','anchor_parent'),('topology','anchor_hsf'),('topology',),()]:
  lookup={}; children=[]; nxt=[]
  for i,m in enumerate(cur):
   k=tuple(m[x] for x in fields)
   if k not in lookup: lookup[k]=len(nxt);nxt.append(m.copy());children.append([])
   children[lookup[k]].append(i)
  code=np.empty(len(cur),dtype=np.int32)
  for i,v in enumerate(children):code[v]=i
  outcomp=np.array([np.unique(comp[v])[0] if len(np.unique(comp[v]))==1 else -1 for v in children],dtype=np.int32)
  levels.append((code,len(nxt),outcomp));cur=nxt;comp=outcomp
 return levels
def aggregate(values,counts,group,levels):
 w=counts[:,group]; v=values
 for code,n,outgroup in levels:
  order=np.argsort(code,kind='stable'); c=code[order]; starts=np.r_[0,np.flatnonzero(c[1:]!=c[:-1])+1]
  den=np.add.reduceat(w[:,order],starts,axis=1); num=np.add.reduceat(v[:,order]*w[:,order],starts,axis=1)
  v=np.divide(num,den,out=np.zeros_like(num),where=den>0);w=(den>0).astype(float)
  same=outgroup>=0
  if same.any():w[:,same]=counts[:,outgroup[same]]
 return v[:,0]
def q(x):return [float(v) for v in np.quantile(x,[.025,.975],method='linear')]
def write_individual(path,rows):
 fields=['bundle_id','anchor','positive','negative','score_anchor_positive','score_anchor_negative','score_difference','F_credit','topology','anchor_hsf','negative_hsf','anchor_parent','positive_parent','negative_parent','point_group','uncertainty_group']
 Path(path).parent.mkdir(parents=True,exist_ok=True)
 with gzip.open(path,'wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rows:w.writerow({k:r[k] for k in fields})
def main():
 p=argparse.ArgumentParser();p.add_argument('--relations',type=Path,required=True);p.add_argument('--objects',type=Path,required=True);p.add_argument('--point-membership',type=Path,required=True);p.add_argument('--uncertainty-membership',type=Path,required=True);p.add_argument('--leaf-scores',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--individual-scores',type=Path,required=True);a=p.parse_args()
 rel=rd(a.relations); obj={r['domain_id']:r for r in rd(a.objects)}; point={r['bundle_id']:r['group_id'] for r in rd(a.point_membership)}; unc={r['bundle_id']:r['group_id'] for r in rd(a.uncertainty_membership)}; leaf={r['bundle_id']:r for r in rd(a.leaf_scores)}
 ids=[r['bundle_id'] for r in rel]
 if len(ids)!=len(set(ids)) or set(ids)!=set(point) or set(ids)!=set(unc) or set(ids)!=set(leaf):raise SystemExit('formal bundle identity mismatch')
 if point!=unc:raise SystemExit('frozen L point and uncertainty memberships disagree')
 rows=[]
 for r in rel:
  l=leaf[r['bundle_id']]; ap=float(l['score_anchor_positive']);an=float(l['score_anchor_negative'])
  if not math.isfinite(ap) or not math.isfinite(an):raise SystemExit('nonfinite leaf score')
  d=ap-an; credit=1.0 if d>FORMAL_DELTA_NUM else -1.0 if d < -FORMAL_DELTA_NUM else 0.0
  for k in ('anchor','positive'): 
   if r[k] not in obj:raise SystemExit('object missing')
  rows.append({'bundle_id':r['bundle_id'],'anchor':r['anchor'],'positive':r['positive'],'negative':r['negative'],'topology':r['topology'],'anchor_hsf':r['anchor_hsf'],'negative_hsf':r['negative_hsf'],'anchor_parent':obj[r['anchor']]['parent_proxy_pdb_chain'],'positive_parent':obj[r['positive']]['parent_proxy_pdb_chain'],'negative_parent':r['negative_parent'],'score_anchor_positive':ap,'score_anchor_negative':an,'score_difference':d,'F_credit':credit,'point_group':point[r['bundle_id']],'uncertainty_group':unc[r['bundle_id']],'value':credit})
 av=anchor_values(rows); anchors=sorted(av); group_names=sorted(set(point.values())); gi={x:i for i,x in enumerate(group_names)}; groups=[]; meta=[]
 for x in anchors:
  gs={point[r['bundle_id']] for r in rows if r['anchor']==x}
  if len(gs)!=1:raise SystemExit('anchor maps to multiple frozen bootstrap groups')
  z=av[x];groups.append(gi[gs.pop()]);meta.append({'topology':z['topology'],'anchor_hsf':z['anchor_hsf'],'anchor_parent':z['anchor_parent']})
 levels=hierarchy(meta,groups); vals=np.array([av[x]['value'] for x in anchors],dtype=float); n=len(group_names); one=np.ones((1,n)); F=float(aggregate(vals[None,:],one,np.array(groups),levels)[0])
 rng=np.random.Generator(np.random.PCG64(FORMAL_SEED)); boot=[]
 for start in range(0,FORMAL_DRAWS,500):
  k=min(500,FORMAL_DRAWS-start); cnt=rng.multinomial(n,np.full(n,1/n),size=k).astype(float); boot.append(aggregate(np.broadcast_to(vals,(k,len(vals))),cnt,np.array(groups),levels))
 boot=np.concatenate(boot); ci=q(boot); out={'artifact_type':'PS02B_FORMAL_SCORER_RESULT','status':'PASS','formal_result':True,'formal_interface':'F_ONLY','score_transform_version':'PS02B_v1.0.2_F_ONLY','F':F,'n_bundles':len(ids),'n_anchors':len(anchors),'n_bootstrap_units':n,'coverage':{'planned_bundles':len(rel),'scored_bundles':len(rows),'fraction':1.0},'delta_num':FORMAL_DELTA_NUM,'bootstrap':{'rng':'PCG64','seed':FORMAL_SEED,'draws':FORMAL_DRAWS,'F_95_CI':ci},'inputs':{'relations_sha256':sha(a.relations),'objects_sha256':sha(a.objects),'point_membership_sha256':sha(a.point_membership),'uncertainty_membership_sha256':sha(a.uncertainty_membership),'leaf_scores_sha256':sha(a.leaf_scores)}}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');write_individual(a.individual_scores,rows)
if __name__=='__main__':main()
