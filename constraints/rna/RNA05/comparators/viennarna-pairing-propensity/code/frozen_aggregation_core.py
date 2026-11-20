#!/usr/bin/env python3
"""Formal F-only sealed-ledger aggregation core."""
import hashlib,json,math
from pathlib import Path
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)
def validate_members(ms,fields):
 ids=set(); member_ids=set(); maps=[{} for _ in fields[1:]]
 for m in ms:
  if m['leaf_id'] in ids:raise ValueError('duplicate leaf identity')
  ids.add(m['leaf_id'])
  key=canon(m)
  if key in member_ids:raise ValueError('duplicate canonical member identity')
  member_ids.add(key)
  for i,(child,parent) in enumerate(zip(fields,fields[1:])):
   if child not in m or parent not in m:raise ValueError('incomplete frozen hierarchy membership')
   if m[child] in maps[i] and maps[i][m[child]]!=m[parent]:raise ValueError('child maps to multiple parent identities')
   maps[i][m[child]]=m[parent]
def safe_cosine(a,b):
 a=[float(x) for x in a];b=[float(x) for x in b]
 if not a or not b:raise ValueError('empty vector')
 if len(a)!=len(b):raise ValueError('vector dimension mismatch')
 if not all(math.isfinite(x) for x in a+b):raise ValueError('nonfinite vector')
 na=math.sqrt(math.fsum(x*x for x in a));nb=math.sqrt(math.fsum(x*x for x in b))
 return 0. if na==0 or nb==0 else math.fsum(x*y for x,y in zip(a,b))/(na*nb)
def derive_credit(r,mode,tau,representation_dimension=None):
 if mode=='sealed_credit':
  v=float(r['F_credit'])
  if not math.isfinite(v) or v not in (-1.,0.,1.):raise ValueError('invalid F-only credit')
  return v
 if mode=='safe_cosine':
  if not isinstance(representation_dimension,int) or representation_dimension<=0:raise ValueError('invalid frozen representation dimension')
  vectors=[r['anchor_vector'],r['positive_vector'],r['negative_vector']]
  if any(len(v)!=representation_dimension for v in vectors):raise ValueError('representation dimension differs from frozen ModelInterface')
  d=safe_cosine(vectors[0],vectors[1])-safe_cosine(vectors[0],vectors[2])
  return 1. if d>tau else (0. if d<-tau else .5)
 raise ValueError('unknown measurement mode')
def aggregate(rows,fields,graph_levels,population_digest):
 ns=[{'id':'EVALUATOR_LEAF_VALUES:'+r['member']['leaf_id'],'v':r['F_credit'],'members':[r['member']]} for r in rows]
 trace=[{'node_id':n['id'],'level':'EVALUATOR_LEAF_VALUES','group_keys':{'leaf_id':n['id'].split(':',1)[1]},'child_ids':[],'child_values':[],'weights':[],'value':n['v'],'population_digest':population_digest} for n in ns]
 for i,f in enumerate(fields):
  level=graph_levels[i];g={}
  for n in ns:
   vals={m[f] for m in n['members']}
   if len(vals)!=1:raise ValueError('aggregate child maps to multiple '+f)
   g.setdefault(next(iter(vals)),[]).append(n)
  ns=[{'id':level+':'+k,'key':k,'v':math.fsum(n['v'] for n in x)/len(x),'members':[m for n in x for m in n['members']],'children':x} for k,x in sorted(g.items())]
  trace.extend({'node_id':n['id'],'level':level,'group_keys':{f:n['key']},'child_ids':[x['id'] for x in n['children']],'child_values':[x['v'] for x in n['children']],'weights':[1/len(n['children'])]*len(n['children']),'value':n['v'],'population_digest':population_digest} for n in ns)
 overall=math.fsum(n['v'] for n in ns)/len(ns)
 trace.append({'node_id':'OVERALL:overall','level':graph_levels[-1],'group_keys':{},'child_ids':[n['id'] for n in ns],'child_values':[n['v'] for n in ns],'weights':[1/len(ns)]*len(ns),'value':overall,'population_digest':population_digest})
 return overall,trace
