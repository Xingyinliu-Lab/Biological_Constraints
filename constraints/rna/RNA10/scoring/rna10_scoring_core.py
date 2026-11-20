#!/usr/bin/env python3
"""C-RNA-10 v0.1.5 fail-closed formal scorer with final-RefSet and source-materialization lineage binding."""
from __future__ import annotations
import hashlib, json, math, random
from collections import defaultdict
from pathlib import Path
def chash(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def avg(v):
    if not v:raise ValueError('empty planned population')
    return sum(v)/len(v)
def load(p):
    x=json.loads(Path(p).read_text())
    if not isinstance(x,dict):raise ValueError('measurement ledger must be an object')
    return x
def validate_measurement_header(plan,ledger,binding,runplan_sha):
    if ledger.get('artifact_type')!='C_RNA10_MEASUREMENT_LEDGER':raise ValueError('wrong measurement ledger artifact_type')
    if not isinstance(plan.get('model_interface_id'),str) or not plan['model_interface_id']:raise ValueError('missing RunPlan model identity')
    if not isinstance(ledger.get('model_id'),str) or not ledger['model_id']:raise ValueError('missing measurement model identity')
    if ledger.get('model_id')!=plan.get('model_interface_id'):raise ValueError('measurement ledger model identity mismatch')
    if not isinstance(ledger.get('measurement_artifact_release_id'),str) or not ledger['measurement_artifact_release_id']:raise ValueError('missing measurement artifact release identity')
    for k,v in (('runplan_sha256',runplan_sha),('population_identity_sha256',binding['population_identity_sha256']),('numeric_method_sha256',binding['numeric_method_sha256']),('split_release_sha256',binding['split_release_sha256'])):
        if ledger.get(k)!=v:raise ValueError('measurement lineage mismatch: '+k)
    if not isinstance(ledger.get('representation_dimension'),int) or ledger['representation_dimension']<=0:raise ValueError('invalid representation dimension')
    if ledger.get('delta_definition')!='z_mut-z_WT':raise ValueError('measurement delta definition mismatch')
    if ledger.get('label_encoding')!={'PRESERVING':0,'DISRUPTING':1}:raise ValueError('measurement label encoding mismatch')
    records=ledger.get('records')
    if not isinstance(records,list):raise ValueError('measurement ledger records missing')
    return records
def validate(plan,b):
    if plan.get('artifact_type')!='C_RNA10_RUNPLAN':raise ValueError('wrong RunPlan artifact_type')
    for k in ('constraint_id','assay_refset_release_id','assay_refset_release_sha256','source_materialization_release_id','source_materialization_manifest_sha256','population_identity_sha256','numeric_method_release_id','numeric_method_sha256','aggregation_contract_sha256','split_release_id','split_release_sha256'):
        if plan.get(k)!=b.get(k):raise ValueError('RunPlan authority mismatch: '+k)
    members=plan.get('canonical_members')
    if not isinstance(members,list) or not members:raise ValueError('empty planned population')
    lids=[m.get('leaf_id') for m in members]
    if any(not x for x in lids) or len(lids)!=len(set(lids)):raise ValueError('duplicate or missing planned leaf_id')
    expected=chash(sorted(lids))
    if plan.get('planned_leaf_set_hash')!=expected:raise ValueError('planned leaf set hash mismatch')
    tau=plan.get('tie_tolerance');draws=plan.get('bootstrap_draws')
    for k in ('tie_tolerance','bootstrap_draws','bootstrap_seed'):
        if plan.get(k)!=b.get(k):raise ValueError('RunPlan Numeric Method constant mismatch: '+k)
    if not isinstance(tau,(int,float)) or not math.isfinite(float(tau)) or float(tau)<0:raise ValueError('invalid tie tolerance')
    if not isinstance(draws,int) or draws<=0:raise ValueError('invalid bootstrap draws')
    if not isinstance(plan.get('bootstrap_seed'),int):raise ValueError('invalid bootstrap seed')
    mapping={}
    for m in members:
        s,u=m.get('source_dependency_unit_id'),m.get('uncertainty_dependency_cluster_id')
        if not s or not u:raise ValueError('missing hierarchy identifier')
        if s in mapping and mapping[s]!=u:raise ValueError('source unit maps to multiple uncertainty clusters')
        mapping[s]=u
    if len(set(mapping.values()))!=len(mapping):raise ValueError('V1 uncertainty membership is not one-to-one')
    hierarchy=sorted((m['leaf_id'],m['matched_block_id'],m['source_dependency_unit_id'],m['uncertainty_dependency_cluster_id']) for m in members)
    if b.get('planned_hierarchy_hash')!=chash(hierarchy):raise ValueError('planned hierarchy hash mismatch')
    return {m['leaf_id']:m for m in members}
def evaluate(plan,ledger,binding,include_trace=False):
    members=validate(plan,binding);runplan_sha=hashlib.sha256(json.dumps(plan,sort_keys=True).encode()).hexdigest();records=validate_measurement_header(plan,ledger,binding,runplan_sha);seen=set();by_block=defaultdict(list);tau=float(plan['tie_tolerance'])
    for r in records:
        lid=r.get('leaf_id')
        if lid in seen:raise ValueError('duplicate measurement leaf')
        seen.add(lid)
        if lid not in members:raise ValueError('unplanned measurement leaf')
        m=members[lid]
        for k in ('matched_block_id','source_dependency_unit_id','uncertainty_dependency_cluster_id'):
            if r.get(k)!=m.get(k):raise ValueError('measurement hierarchy mismatch: '+k)
        try:hp,hn=float(r['h_positive']),float(r['h_negative'])
        except (TypeError,ValueError):raise ValueError('invalid measurement score')
        if not math.isfinite(hp) or not math.isfinite(hn):raise ValueError('non-finite measurement score')
        F_credit=1.0 if hp-hn>tau else -1.0 if hn-hp>tau else 0.0
        by_block[m['matched_block_id']].append((m,F_credit))
    if seen!=set(members):raise ValueError('measurement ledger does not close planned population')
    block_score={};source_blocks=defaultdict(list);source_u={}
    for block,leaves in by_block.items():
        ss={x[0]['source_dependency_unit_id'] for x in leaves};uu={x[0]['uncertainty_dependency_cluster_id'] for x in leaves}
        if len(ss)!=1 or len(uu)!=1:raise ValueError('child maps to multiple parents')
        s,u=next(iter(ss)),next(iter(uu));source_u[s]=u
        block_score[block]=avg([x[1] for x in leaves]);source_blocks[s].append(block_score[block])
    source_F={s:avg(x) for s,x in source_blocks.items()};point_F=avg(list(source_F.values()))
    u_scores={source_u[s]:v for s,v in source_F.items()}
    rng=random.Random(plan['bootstrap_seed']);ids=sorted(u_scores);draws=sorted(avg([u_scores[rng.choice(ids)] for _ in ids]) for _ in range(plan['bootstrap_draws']))
    lo,hi=draws[int(.025*(len(draws)-1))],draws[int(.975*(len(draws)-1))]
    trace={'leaf_credits':[{'leaf_id':m['leaf_id'],'F_credit':c,'matched_block_id':m['matched_block_id'],'source_dependency_unit_id':m['source_dependency_unit_id'],'uncertainty_dependency_cluster_id':m['uncertainty_dependency_cluster_id']} for b in sorted(by_block) for m,c in by_block[b]],'block_F':block_score,'source_F':source_F,'point_F':point_F,'F_bootstrap':draws}
    result={'F':point_F,'F_95_CI':[lo,hi],'formal_interface':'F_ONLY','score_transform_version':'RNA10_F_ONLY_v1','blocks':len(block_score),'source_dependency_units':len(source_F),'uncertainty_dependency_clusters':len(u_scores),'aggregation_trace_sha256':chash(trace)}
    if include_trace:result['_aggregation_trace']=trace
    return result
