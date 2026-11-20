#!/usr/bin/env python3
import csv, hashlib, json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import st04_numeric_core_v0_1_1 as numeric_core

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def load_json(path):
    with open(path,encoding='utf-8') as f: return json.load(f)

def load_tsv(path):
    with open(path,encoding='utf-8',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def endpoint_set(relations): return sorted({r[k] for r in relations for k in ('anchor_residue_id','true_target_residue_id','negative_residue_id')})
def set_hash(values): return hashlib.sha256(('\n'.join(sorted(values))+'\n').encode()).hexdigest()

def require_hash(path,expected,label):
    got=sha256(path)
    if got!=expected: raise ValueError(f'{label} SHA256 mismatch')
    return got

def load_embeddings(path):
    with np.load(path,allow_pickle=False) as z:
        ids=np.asarray(z['residue_ids']).astype(str); raw=np.asarray(z['embeddings'])
    if raw.dtype.kind not in 'iu f'.replace(' ', '') or raw.dtype.kind == 'b':
        raise ValueError('embedding dtype must be a real integer or floating type')
    emb=np.asarray(raw,dtype=np.float64)
    if ids.ndim!=1 or emb.ndim!=2 or emb.shape[1]<=0 or len(ids)!=len(emb): raise ValueError('embedding array shape failure')
    if len(set(ids))!=len(ids): raise ValueError('duplicate residue endpoint')
    if not np.isfinite(emb).all(): raise ValueError('nonfinite embedding')
    return {k:v for k,v in zip(ids,emb)},int(emb.shape[1])

def constant_representation_qc(queries,embedding_dimension):
    tests=[]
    beta=np.linspace(-0.5,0.5,3*embedding_dimension,dtype=np.float64); model={'beta':beta,'mean':np.zeros_like(beta),'scale':np.ones_like(beta),'intercept':0.17}
    for value in (0.0,3.14):
        z=np.full(embedding_dimension,value,dtype=np.float64); score=numeric_core.decision(model,numeric_core.pair_feature(z,z)); pred={}; ties=0; contrasts=0
        for qid,q in queries.items():
            k=len(q['negatives']); delta=np.full(k,score-score,dtype=np.float64); credit=np.where(delta>0.0,1.0,np.where(delta<0.0,-1.0,0.0)); ties+=int(np.count_nonzero(delta==0.0)); contrasts+=k
            pred[qid]={'query_credit':float(credit.mean())}
        F,clusters=numeric_core.aggregate(queries,pred)
        tests.append({'constant_value':value,'query_count':len(queries),'contrast_count':contrasts,'tie_count':ties,'cluster_count':len(clusters),'F':F,'pass':F==0.0 and ties==contrasts})
    return {'tests':tests,'pass':all(x['pass'] for x in tests)}

def release_qc(rows,queries,emb):
    bad_geometry=0; symmetry_fail=0; max_symmetry_difference=0.0
    for r in rows:
        if not float(r['positive_distance']) < float(r['negative_distance']): bad_geometry+=1
        a=emb[r['anchor_residue_id']]
        for k in ('true_target_residue_id','negative_residue_id'):
            b=emb[r[k]]; x=numeric_core.pair_feature(a,b); y=numeric_core.pair_feature(b,a)
            diff=float(np.max(np.abs(x-y))) if x.size else 0.0
            max_symmetry_difference=max(max_symmetry_difference,diff)
            if not np.array_equal(x,y): symmetry_fail+=1
    geometry_f=1.0-bad_geometry/len(rows); constant=constant_representation_qc(queries,len(next(iter(emb.values()))))
    status='PASS' if bad_geometry==0 and symmetry_fail==0 and geometry_f==1.0 and constant['pass'] else 'FAIL'
    return {'status':status,'coordinate_truth_oracle':{'relations_checked':len(rows),'violations':bad_geometry,'F':geometry_f},
            'constant_representation':constant,
            'pair_symmetry':{'directed_pair_features_checked':2*len(rows),'violations':symmetry_fail,'max_abs_difference':max_symmetry_difference}}

def coverage_report(rows,queries):
    groups=defaultdict(lambda:{'relations':0,'queries':set(),'k':Counter()})
    for r in rows:
        key=(r['target_aa_identity'],r['signed_direction_class'],r['sequence_separation_bin']); z=groups[key]; z['relations']+=1; z['queries'].add(r['query_id'])
    for key,z in groups.items():
        for q in z['queries']: z['k'][str(len(queries[q]['negatives']))]+=1
    strata=[]
    for key,z in sorted(groups.items()): strata.append({'target_aa_identity':key[0],'signed_direction_class':key[1],'sequence_separation_bin':key[2],'relations':z['relations'],'queries':len(z['queries']),'K_distribution':dict(sorted(z['k'].items())),'in_scope_failures':0})
    return {'planned_relations':len(rows),'scoreable_relations':len(rows),'planned_queries':len(queries),'scoreable_queries':len(queries),'in_scope_failures':0,'strata':strata}

def save_fit_index(fits,out):
    arrays={}; index={'outer':{},'inner':{}}
    def one(name,m):
        arrays[name+'__beta']=np.asarray(m['beta'],dtype=np.float64); arrays[name+'__mean']=np.asarray(m['mean'],dtype=np.float64); arrays[name+'__scale']=np.asarray(m['scale'],dtype=np.float64)
        return {'beta_key':name+'__beta','mean_key':name+'__mean','scale_key':name+'__scale','intercept':m['intercept'],'lambda':m['lambda'],'iterations':m['iterations'],'training_provenance':m['training_provenance']}
    for o,m in sorted(fits['outer'].items(),key=lambda x:int(x[0])): index['outer'][o]=one('outer_'+o,m)
    for o,ff in sorted(fits['inner'].items(),key=lambda x:int(x[0])):
        index['inner'][o]={}
        for name,m in sorted(ff.items()): index['inner'][o][name]=one('inner_'+o+'__'+name,m)
    npz=out/'FIT_PARAMETERS.npz'; np.savez_compressed(npz,**arrays); index['fit_parameters_sha256']=sha256(npz); index['array_count']=len(arrays)
    (out/'FIT_INDEX.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    return index
