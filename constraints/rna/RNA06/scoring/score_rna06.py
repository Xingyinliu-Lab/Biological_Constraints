from __future__ import annotations
import argparse, csv, hashlib, json, re
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

LAMBDA_GRID=(1e-4,1e-3,1e-2,1e-1,1.,10.,100.)

INTEGER = re.compile(r"0|[1-9][0-9]*")

def strict_nonnegative_integer(value, field):
    text=str(value)
    if not INTEGER.fullmatch(text):
        raise ValueError(f'{field} must be a non-negative integer: {value!r}')
    return int(text)

def real_matrix(value, label):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf':
        raise ValueError(f'{label} must contain real numeric values')
    matrix=np.asarray(raw,dtype=np.float32)
    if matrix.ndim!=2 or matrix.shape[1]<=0 or not np.isfinite(matrix).all():
        raise ValueError(f'invalid {label}')
    return matrix

def load_vectors(path):
    if str(path).lower().endswith('.npz'):
        d=np.load(path,allow_pickle=False)
        if set(d.files)=={'sequence_sha256','position0','vectors'}:
            sequence_ids=np.asarray(d['sequence_sha256'])
            positions=np.asarray(d['position0'])
            if sequence_ids.ndim!=1 or sequence_ids.dtype.kind not in 'SU': raise ValueError('sequence_sha256 must be a one-dimensional string array')
            if positions.ndim!=1 or positions.dtype.kind not in 'iu': raise ValueError('position0 must be a one-dimensional integer array')
            if np.issubdtype(positions.dtype,np.signedinteger) and np.any(positions<0): raise ValueError('position0 must be non-negative')
            keys=[(str(s),int(p)) for s,p in zip(sequence_ids,positions)]; a=real_matrix(d['vectors'],'position-vector NPZ')
            if len(keys)!=len(a) or len(set(keys))!=len(keys): raise ValueError('invalid position-vector NPZ')
            return {k:v for k,v in zip(keys,a)},'position'
        matrices={}
        for key in d.files:
            value=real_matrix(d[key],f'sequence matrix: {key}')
            matrices[str(key)]=value
        if not matrices: raise ValueError('empty sequence-keyed NPZ')
        return matrices,'matrix_npz'
    import torch
    loaded=torch.load(path,map_location='cpu',weights_only=True)
    if not isinstance(loaded,dict) or not loaded: raise ValueError('PyTorch compatibility input must be a non-empty mapping')
    matrices={}
    for key,value in loaded.items():
        if not isinstance(value,torch.Tensor) or value.dtype in (torch.bool,torch.complex64,torch.complex128):
            raise ValueError(f'invalid sequence matrix: {key}')
        matrices[str(key)]=real_matrix(value.detach().cpu().numpy(),f'sequence matrix: {key}')
    return matrices,'matrix'

def read(path):
    with open(path, newline='', encoding='utf-8') as h:
        return list(csv.DictReader(h, delimiter='\t'))

def read_fasta(path):
    records={}; ident=None; chunks=[]
    for raw in Path(path).read_text(encoding='utf-8').splitlines():
        if raw.startswith('>'):
            if ident is not None:
                if ident in records: raise ValueError(f'duplicate FASTA record: {ident}')
                records[ident]=''.join(chunks).upper()
            ident=raw[1:].split()[0]; chunks=[]
        else:
            chunks.append(raw.strip())
    if ident is not None:
        if ident in records: raise ValueError(f'duplicate FASTA record: {ident}')
        records[ident]=''.join(chunks).upper()
    if not records: raise ValueError('empty RNA FASTA')
    return records

def separation_bin(distance):
    for upper,label in ((4,'01_04'),(8,'05_08'),(16,'09_16'),(32,'17_32'),(64,'33_64'),(128,'65_128')):
        if distance<=upper: return label
    return '129_plus'

def validate_refset(fasta_path, object_path, directed_path, candidate_path, edge_path, pair_state_path, manifest_path, objects, directed, candidates, edges, pair_states):
    object_by_id={}
    for row in objects:
        record=row['record_id']; sequence_id=row['sequence_sha256']
        if record in object_by_id: raise ValueError('duplicate RNA object identity')
        strict_nonnegative_integer(row['length_nt'],'length_nt')
        if not sequence_id: raise ValueError(f'missing sequence identity: {record}')
        object_by_id[record]=row
    sequence_by_record={record:row['sequence_sha256'] for record,row in object_by_id.items()}
    edge_by_id={}
    for edge in edges:
        edge_id=edge['edge_bundle_id']; record=edge['record_id']
        if edge_id in edge_by_id or record not in object_by_id: raise ValueError(f'invalid reciprocal edge: {edge_id}')
        i=strict_nonnegative_integer(edge['position_i0'],'position_i0'); j=strict_nonnegative_integer(edge['position_j0'],'position_j0')
        if i==j: raise ValueError(f'invalid reciprocal edge coordinates: {edge_id}')
        if strict_nonnegative_integer(edge['reciprocal_directed_challenges'],'reciprocal_directed_challenges')!=2 or edge['complete_execution_bundle'].lower()!='true': raise ValueError(f'incomplete reciprocal edge: {edge_id}')
        edge_by_id[edge_id]=edge
    endpoint_by_id={}; endpoints_by_edge=defaultdict(list)
    for endpoint in directed:
        endpoint_id=endpoint['directed_endpoint_id']; edge_id=endpoint['edge_bundle_id']; record=endpoint['record_id']
        if endpoint_id in endpoint_by_id or edge_id not in edge_by_id or record not in object_by_id: raise ValueError(f'invalid directed endpoint: {endpoint_id}')
        anchor=strict_nonnegative_integer(endpoint['anchor_position0'],'anchor_position0'); partner=strict_nonnegative_integer(endpoint['true_partner_position0'],'true_partner_position0')
        if anchor==partner: raise ValueError(f'invalid directed endpoint coordinates: {endpoint_id}')
        edge=edge_by_id[edge_id]
        if record!=edge['record_id'] or {anchor,partner}!={strict_nonnegative_integer(edge['position_i0'],'position_i0'),strict_nonnegative_integer(edge['position_j0'],'position_j0')}: raise ValueError(f'endpoint/edge mismatch: {endpoint_id}')
        for field in ('outer_dependency_cluster_id','parent_id','family_id','primary_structure_map_id','primary_condition_id'):
            if endpoint[field]!=edge[field]: raise ValueError(f'endpoint/edge {field} mismatch: {endpoint_id}')
        expected_role='I_TO_J' if anchor==strict_nonnegative_integer(edge['position_i0'],'position_i0') else 'J_TO_I'
        if endpoint['role']!=expected_role or int(endpoint['signed_direction'])!=(1 if partner>anchor else -1): raise ValueError(f'invalid endpoint direction: {endpoint_id}')
        distance=abs(partner-anchor)
        if int(endpoint['separation'])!=distance or endpoint['separation_bin']!=separation_bin(distance): raise ValueError(f'invalid endpoint separation: {endpoint_id}')
        endpoint_by_id[endpoint_id]=endpoint; endpoints_by_edge[edge_id].append(endpoint)
    if set(endpoints_by_edge)!=set(edge_by_id): raise ValueError('directed endpoints do not cover reciprocal edges')
    for edge_id,rows in endpoints_by_edge.items():
        if len(rows)!=2:
            raise ValueError(f'reciprocal edge does not have two endpoints: {edge_id}')
        pairs={(int(row['anchor_position0']),int(row['true_partner_position0'])) for row in rows}
        edge=edge_by_id[edge_id]; i=int(edge['position_i0']); j=int(edge['position_j0'])
        if pairs!={(i,j),(j,i)}: raise ValueError(f'non-reciprocal endpoints: {edge_id}')
    candidates_by_endpoint=defaultdict(list); candidate_ids=set()
    for candidate in candidates:
        candidate_id=candidate['candidate_id']; endpoint_id=candidate['directed_endpoint_id']
        if candidate_id in candidate_ids or endpoint_id not in endpoint_by_id: raise ValueError(f'invalid candidate: {candidate_id}')
        candidate_ids.add(candidate_id); endpoint=endpoint_by_id[endpoint_id]; record=endpoint['record_id']
        for field in ('edge_bundle_id','record_id','anchor_position0','true_partner_position0','outer_dependency_cluster_id','selection_contract_id'):
            if candidate[field]!=endpoint[field]: raise ValueError(f'candidate/endpoint {field} mismatch: {candidate_id}')
        partner=strict_nonnegative_integer(candidate['paired_nonpartner_position0'],'paired_nonpartner_position0')
        if partner in {int(endpoint['anchor_position0']),int(endpoint['true_partner_position0'])}: raise ValueError(f'invalid nonpartner coordinate: {candidate_id}')
        strict_nonnegative_integer(candidate['candidate_degree'],'candidate_degree')
        anchor=int(endpoint['anchor_position0']); true_partner=int(endpoint['true_partner_position0'])
        if int(candidate['signed_direction'])!=(1 if partner>anchor else -1): raise ValueError(f'candidate direction mismatch: {candidate_id}')
        if int(candidate['true_separation'])!=abs(true_partner-anchor) or int(candidate['candidate_separation'])!=abs(partner-anchor): raise ValueError(f'candidate separation mismatch: {candidate_id}')
        if candidate['separation_bin']!=endpoint['separation_bin']: raise ValueError(f'candidate separation bin mismatch: {candidate_id}')
        candidates_by_endpoint[endpoint_id].append(candidate)
    if set(candidates_by_endpoint)!=set(endpoint_by_id): raise ValueError('candidate ledger does not cover directed endpoints')
    for endpoint_id,rows in candidates_by_endpoint.items():
        endpoint=endpoint_by_id[endpoint_id]; selected=strict_nonnegative_integer(endpoint['selected_K_q'],'selected_K_q')
        if len(rows)!=selected or any(strict_nonnegative_integer(row['selected_K_q'],'selected_K_q')!=selected for row in rows): raise ValueError(f'candidate count mismatch: {endpoint_id}')
        if len({row['paired_nonpartner_position0'] for row in rows})!=selected or {strict_nonnegative_integer(row['selection_rank'],'selection_rank') for row in rows}!=set(range(1,selected+1)): raise ValueError(f'candidate selection identity mismatch: {endpoint_id}')
    manifest=json.loads(Path(manifest_path).read_text(encoding='utf-8'))
    expected={'directed_endpoint_challenges':len(directed),'paired_nonpartner_candidates_after_cap':len(candidates),'edge_bundles':len(edges),'parents':len({r['parent_id'] for r in directed}),'outer_dependency_clusters':len({r['outer_dependency_cluster_id'] for r in directed})}
    for field,value in expected.items():
        if manifest.get(field)!=value: raise ValueError(f'population manifest {field} mismatch')
    return sequence_by_record, {'records':len(objects),'unique_sequences':len(set(sequence_by_record.values())),**expected,'fasta_sha256':sha256(fasta_path),'objects_sha256':sha256(object_path),'directed_sha256':sha256(directed_path),'candidates_sha256':sha256(candidate_path),'edges_sha256':sha256(edge_path),'pair_states_sha256':sha256(pair_state_path),'manifest_sha256':sha256(manifest_path)}

def folds(groups, n=4):
    ordered=sorted(set(groups), key=lambda x: hashlib.sha256(x.encode()).hexdigest())
    return [{x for j,x in enumerate(ordered) if j % n == i} for i in range(n)]

def build_samples(directed, candidates, sequence_by_record, vectors, vector_mode):
    by_endpoint=defaultdict(list)
    for r in candidates: by_endpoint[r['directed_endpoint_id']].append(r)
    samples=[]
    for e in directed:
        endpoint=e['directed_endpoint_id']; seq=sequence_by_record[e['record_id']]
        anchor=int(e['anchor_position0']); true=int(e['true_partner_position0'])
        if endpoint not in by_endpoint: raise ValueError(f'missing candidates: {endpoint}')
        def feature(partner):
            if vector_mode=='position':
                if (seq,anchor) not in vectors or (seq,partner) not in vectors: raise ValueError('missing nucleotide position vector')
                z_a=np.asarray(vectors[(seq,anchor)],dtype=np.float32); z_b=np.asarray(vectors[(seq,partner)],dtype=np.float32)
            else:
                if seq not in vectors: raise ValueError(f'missing vector sequence: {seq}')
                if anchor<0 or partner<0 or anchor>=vectors[seq].shape[0] or partner>=vectors[seq].shape[0]: raise ValueError('invalid nucleotide coordinate mapping')
                z_a=np.asarray(vectors[seq][anchor],dtype=np.float32)
                z_b=np.asarray(vectors[seq][partner],dtype=np.float32)
            if z_a.ndim != 1 or z_b.ndim != 1 or z_a.shape != z_b.shape: raise ValueError('invalid nucleotide vectors')
            if not (np.isfinite(z_a).all() and np.isfinite(z_b).all()): raise ValueError('nonfinite nucleotide vectors')
            return np.concatenate((z_a+z_b,z_a*z_b))
        common={k:e[k] for k in ('outer_dependency_cluster_id','parent_id','edge_bundle_id','directed_endpoint_id')}
        samples.append({**common,'sample_id':endpoint+':true','candidate_id':'','anchor_position0':anchor,'partner_position0':true,'class':'true_partner','label':1.,'feature':feature(true)})
        for c in by_endpoint[endpoint]:
            partner=int(c['paired_nonpartner_position0'])
            samples.append({**common,'sample_id':c['candidate_id'],'candidate_id':c['candidate_id'],'anchor_position0':anchor,'partner_position0':partner,'class':'paired_nonpartner','label':0.,'feature':feature(partner)})
    dim={s['feature'].shape[0] for s in samples}
    if len(dim)!=1: raise ValueError(f'feature dimensions differ: {dim}')
    return samples, np.stack([s.pop('feature') for s in samples]).astype(np.float32,copy=False)

def hierarchy_weights(rows):
    # cluster -> parent -> edge -> endpoint -> class -> row, equal at every level
    tree=defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(list)))))
    for i,r in enumerate(rows):
        tree[r['outer_dependency_cluster_id']][r['parent_id']][r['edge_bundle_id']][r['directed_endpoint_id']][r['class']].append(i)
    w=np.zeros(len(rows),float)
    for parents in tree.values():
        for edges in parents.values():
            for endpoints in edges.values():
                for classes in endpoints.values():
                    for ids in classes.values():
                        w[ids]=1/len(tree)/len(parents)/len(edges)/len(endpoints)/len(classes)/len(ids)
    return w/w.mean()

def aggregate(rows, scores):
    by_endpoint=defaultdict(lambda: defaultdict(list)); endpoint_meta={}
    for r,s in zip(rows,scores):
        by_endpoint[r['directed_endpoint_id']][r['class']].append(float(s))
        endpoint_meta[r['directed_endpoint_id']] = (r['edge_bundle_id'],r['parent_id'],r['outer_dependency_cluster_id'])
    endpoint_value={}; endpoint_raw={}
    for endpoint,classes in by_endpoint.items():
        if len(classes['true_partner']) != 1 or not classes['paired_nonpartner']:
            raise ValueError(f'incomplete endpoint: {endpoint}')
        h_true=classes['true_partner'][0]
        credits=[1. if h_true>x else (-1. if h_true<x else 0.) for x in classes['paired_nonpartner']]
        endpoint_value[endpoint]=float(np.mean(credits))
        endpoint_raw[endpoint]={'true_score':float(h_true),'nonpartner_mean_score':float(np.mean(classes['paired_nonpartner'])),'nonpartner_count':len(classes['paired_nonpartner'])}
    by_edge=defaultdict(list); edge_meta={}
    for endpoint,value in endpoint_value.items():
        edge,parent,cluster=endpoint_meta[endpoint]
        by_edge[edge].append(value); edge_meta[edge]=(parent,cluster)
    edge_value={edge:float(np.mean(values)) for edge,values in by_edge.items()}
    parent_values=defaultdict(list); parent_cluster={}
    for edge,values in by_edge.items():
        if len(values)!=2: raise ValueError(f'incomplete reciprocal edge: {edge}')
        parent,cluster=edge_meta[edge]; parent_values[parent].append(float(np.mean(values))); parent_cluster[parent]=cluster
    cluster_values=defaultdict(list)
    for parent,values in parent_values.items(): cluster_values[parent_cluster[parent]].append(float(np.mean(values)))
    parent_score={parent:float(np.mean(values)) for parent,values in parent_values.items()}
    cluster_score={cluster:float(np.mean(values)) for cluster,values in cluster_values.items()}
    return float(np.mean(list(cluster_score.values()))),cluster_score,endpoint_value,endpoint_raw,edge_value,parent_score,parent_cluster

def macro(rows, scores):
    F,clusters,*_=aggregate(rows,scores)
    return F,clusters

def fit_predict(train_idx,test_idx,X,rows,lam):
    train_rows=[rows[i] for i in train_idx]
    w=hierarchy_weights(train_rows); w=w/w.sum()
    x=X[train_idx].astype(float,copy=False)
    mu=(x*w[:,None]).sum(0)
    sd=np.sqrt(((x-mu)**2*w[:,None]).sum(0)); sd[sd==0]=1.
    xx=(x-mu)/sd; xt=(X[test_idx].astype(float,copy=False)-mu)/sd
    y=np.array([rows[i]['label'] for i in train_idx],float)
    def objective(theta):
        z=xx@theta[:-1]+theta[-1]
        p=1/(1+np.exp(-np.clip(z,-40,40)))
        d=(p-y)*w
        loss=float(np.sum(w*(np.logaddexp(0,z)-y*z))+.5*lam*np.dot(theta[:-1],theta[:-1]))
        grad=np.r_[xx.T@d+lam*theta[:-1],d.sum()]
        return loss,grad
    res=minimize(lambda q: objective(q),np.zeros(xx.shape[1]+1),jac=True,method='L-BFGS-B',options={'maxiter':2000})
    if not res.success or not np.isfinite(res.x).all(): raise RuntimeError(res.message)
    return xt@res.x[:-1]+res.x[-1]

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    root=Path(__file__).resolve().parents[1]
    ap.add_argument('--directed',default=str(root/'refset/RNA06_FORMAL_DIRECTED_ENDPOINT_LEDGER.tsv'))
    ap.add_argument('--candidates',default=str(root/'refset/RNA06_FORMAL_PAIRED_NONPARTNER_CANDIDATE_LEDGER.tsv'))
    ap.add_argument('--edges',default=str(root/'refset/RNA06_FORMAL_RECIPROCAL_EDGE_BUNDLE_LEDGER.tsv'))
    ap.add_argument('--pair-states',default=str(root/'refset/RNA06_FROZEN_PAIR_STATE_MAPS.tsv'))
    ap.add_argument('--sequences',default=str(root/'refset/rna_sequences.fasta'))
    ap.add_argument('--population-manifest',default=str(root/'refset/PLANNED_POPULATION_MANIFEST.json'))
    source=ap.add_mutually_exclusive_group(required=False)
    source.add_argument('--objects',help='RNA06 rna_objects.tsv (recommended public interface)')
    source.add_argument('--rna05-endpoints',help='historical compatibility input')
    ap.add_argument('--vectors',required=True); ap.add_argument('--model',required=True)
    ap.add_argument('--output',help='output directory')
    ap.add_argument('--out',help=argparse.SUPPRESS)
    ap.add_argument('--bootstrap',type=int,default=20000,help=argparse.SUPPRESS)
    a=ap.parse_args()
    if a.bootstrap != 20000: raise ValueError('formal RNA06 scoring fixes bootstrap draws at 20000')
    if not (a.output or a.out): raise ValueError('--output is required')
    output=Path(a.output or a.out)
    object_path=a.objects or a.rna05_endpoints or str(root/'refset/rna_objects.tsv')
    directed=read(a.directed); candidates=read(a.candidates); objects=read(object_path); edges=read(a.edges); pair_states=read(a.pair_states)
    sequence_by_record,refset_identity=validate_refset(a.sequences,object_path,a.directed,a.candidates,a.edges,a.pair_states,a.population_manifest,objects,directed,candidates,edges,pair_states)
    vectors,vector_mode=load_vectors(a.vectors)
    if vector_mode=='matrix_npz' and set(vectors)!=set(sequence_by_record.values()):
        raise ValueError('sequence-keyed NPZ must exactly cover the released RNA sequence set')
    rows,X=build_samples(directed,candidates,sequence_by_record,vectors,vector_mode)
    groups=np.array([r['outer_dependency_cluster_id'] for r in rows]); pred=np.full(len(rows),np.nan); selected=[]
    for held in sorted(set(groups)):
        train=np.where(groups != held)[0]; test=np.where(groups == held)[0]; choice=None
        for lam in LAMBDA_GRID:
            inner=[]
            for held_inner in folds(groups[train]):
                inner_train=train[np.array([g not in held_inner for g in groups[train]])]
                inner_test=train[np.array([g in held_inner for g in groups[train]])]
                inner.append(macro([rows[i] for i in inner_test],fit_predict(inner_train,inner_test,X,rows,lam))[0])
            key=(float(np.mean(inner)),-lam)
            if choice is None or key>choice[0]: choice=(key,lam)
        pred[test]=fit_predict(train,test,X,rows,choice[1]); selected.append({'held_out_cluster':held,'lambda':choice[1],'inner_F':choice[0][0]})
    F,cluster_scores,endpoint_scores,endpoint_raw,edge_scores,parent_scores,parent_cluster=aggregate(rows,pred)
    rng=np.random.default_rng(20261003); cluster=np.array(list(cluster_scores.values()))
    boot=rng.choice(cluster,(a.bootstrap,len(cluster)),replace=True).mean(1); lo,hi=np.quantile(boot,[.025,.975])
    output.mkdir(parents=True,exist_ok=True)
    pred_rows=[{**{k:r[k] for k in ('sample_id','candidate_id','directed_endpoint_id','edge_bundle_id','parent_id','outer_dependency_cluster_id','anchor_position0','partner_position0','class','label')},'heldout_score':float(s)} for r,s in zip(rows,pred)]
    write_tsv(output/'HELDOUT_PREDICTIONS.tsv',pred_rows,list(pred_rows[0]))
    true_score={r['directed_endpoint_id']:float(s) for r,s in zip(rows,pred) if r['class']=='true_partner'}
    true_position={r['directed_endpoint_id']:r['partner_position0'] for r in rows if r['class']=='true_partner'}
    candidate_rows=[]
    for r,s in zip(rows,pred):
        if r['class']!='paired_nonpartner': continue
        t=true_score[r['directed_endpoint_id']]; n=float(s)
        candidate_rows.append({'candidate_id':r['candidate_id'],'directed_endpoint_id':r['directed_endpoint_id'],'edge_bundle_id':r['edge_bundle_id'],'parent_id':r['parent_id'],'outer_dependency_cluster_id':r['outer_dependency_cluster_id'],'anchor_position0':r['anchor_position0'],'true_partner_position0':true_position[r['directed_endpoint_id']],'paired_nonpartner_position0':r['partner_position0'],'true_score':t,'nonpartner_score':n,'margin':t-n,'F_credit':1.0 if t>n else (-1.0 if t<n else 0.0)})
    write_tsv(output/'CANDIDATE_SCORES.tsv',candidate_rows,list(candidate_rows[0]))
    endpoint_rows=[]
    endpoint_meta={r['directed_endpoint_id']:r for r in rows}
    for endpoint,value in sorted(endpoint_scores.items()):
        meta=endpoint_meta[endpoint]; raw=endpoint_raw[endpoint]
        endpoint_rows.append({'directed_endpoint_id':endpoint,'edge_bundle_id':meta['edge_bundle_id'],'parent_id':meta['parent_id'],'outer_dependency_cluster_id':meta['outer_dependency_cluster_id'],**raw,'individual_F':value})
    write_tsv(output/'INDIVIDUAL_SCORES.tsv',endpoint_rows,list(endpoint_rows[0]))
    edge_meta={r['edge_bundle_id']:r for r in rows}
    edge_rows=[{'edge_bundle_id':edge,'parent_id':edge_meta[edge]['parent_id'],'outer_dependency_cluster_id':edge_meta[edge]['outer_dependency_cluster_id'],'edge_F':score} for edge,score in sorted(edge_scores.items())]
    write_tsv(output/'EDGE_SCORE_LEDGER.tsv',edge_rows,list(edge_rows[0]))
    parent_rows=[{'parent_id':parent,'outer_dependency_cluster_id':parent_cluster[parent],'parent_F':score} for parent,score in sorted(parent_scores.items())]
    write_tsv(output/'PARENT_SCORE_LEDGER.tsv',parent_rows,list(parent_rows[0]))
    cluster_rows=[{'outer_dependency_cluster_id':cluster_id,'cluster_F':score} for cluster_id,score in sorted(cluster_scores.items())]
    write_tsv(output/'CLUSTER_SCORE_LEDGER.tsv',cluster_rows,list(cluster_rows[0]))
    np.save(output/'BOOTSTRAP_MEANS.npy',boot)
    result={'constraint':'RNA06','release':'RNA06_v0.1.3_F_ONLY','formal_result':True,'formal_interface':'RNA06_F_ONLY_v0.1.3','score_transform_version':'direct_F_credit_v1','model':a.model,'F':F,'F_95_CI':[float(lo),float(hi)],'coverage':{'clusters':len(cluster_scores),'parents':len(parent_scores),'edges':len(edge_scores),'directed_endpoints':len(directed),'candidate_rows':len(candidates),'prediction_rows':len(rows)},'outer_LOO':selected,'bootstrap':{'draws':20000,'seed':20261003,'unit':'outer_dependency_cluster_F','method':'PCG64 percentile linear'},'lambda_grid':list(LAMBDA_GRID),'feature':'sum_plus_hadamard_linear_ridge','inputs':{'refset':refset_identity,'vectors_sha256':sha256(a.vectors)}}
    (output/'FORMAL_RESULT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    files=[p for p in output.iterdir() if p.is_file() and p.name!='OUTPUT_MANIFEST.json']
    manifest={'formal_interface':'RNA06_F_ONLY_v0.1.3','files':[{'path':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)} for p in sorted(files)]}
    (output/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__': main()
