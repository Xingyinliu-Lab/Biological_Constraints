#!/usr/bin/env python3
"""Formal PS04 population gate; rejects any runtime leaf-population shrink."""
import csv,gzip,hashlib

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def _open(p): return gzip.open(p,'rt',newline='') if str(p).endswith('.gz') else open(p,newline='')
def relation_identity(r): return (r['anchor_id'],r['positive_id'],r['negative_id'],r['negative_hsf_id'],r['negative_topology_id'],r['positive_hsf_id'],r['anchor_parent_proxy'],r['anchor_source_unit_id'],r['anchor_hsf_id'],r['anchor_topology_id'],r['anchor_architecture_id'])
def leaf_identity(r): return (r['anchor_id'],r['positive_endpoint_id'],r['negative_endpoint_id'],r['negative_hsf_id'],r['negative_topology_id'],r['positive_hsf_id'],r['anchor_parent_id'],r['anchor_source_unit_id'],r['anchor_hsf_id'],r['anchor_topology_id'],r['anchor_architecture_id'])
def assert_complete_population(formal_relations,leaf_scores):
 n=0
 with _open(formal_relations) as a,_open(leaf_scores) as b:
  ra,rb=csv.DictReader(a,delimiter='\t'),csv.DictReader(b,delimiter='\t')
  for n,(x,y) in enumerate(zip(ra,rb),start=1):
   if relation_identity(x)!=leaf_identity(y): raise ValueError(f'formal leaf mismatch at row {n}')
  if next(ra,None) is not None or next(rb,None) is not None: raise ValueError('formal leaf cardinality mismatch')
 return {'status':'PASS','relation_count':n,'formal_relation_ledger_sha256':sha(formal_relations),'leaf_ledger_sha256':sha(leaf_scores)}
