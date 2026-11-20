#!/usr/bin/env python3
from collections import defaultdict
import hashlib
import numpy as np
from scipy.optimize import minimize

def _hash_lines(values):
    return hashlib.sha256(('\n'.join(map(str,values))+'\n').encode()).hexdigest()

def _hash_f64(values):
    return hashlib.sha256(np.asarray(values,dtype='<f8').tobytes(order='C')).hexdigest()

def pair_feature(a,b):
    a=np.asarray(a,dtype=np.float64); b=np.asarray(b,dtype=np.float64)
    return np.concatenate((a+b,np.abs(a-b),a*b))

def build_queries(relations):
    q={}
    for r in relations:
        k=r['query_id']; base=(r['source_dependency_cluster_id'],r['structure_object_id'],r['anchor_residue_id'],r['true_target_residue_id'])
        if k not in q: q[k]={'cluster':base[0],'structure':base[1],'anchor':base[2],'true':base[3],'negatives':[]}
        if (q[k]['cluster'],q[k]['structure'],q[k]['anchor'],q[k]['true'])!=base: raise ValueError('query endpoint closure failure')
        q[k]['negatives'].append((int(r['within_query_rank']),r['negative_residue_id']))
    for k,v in q.items():
        v['negatives'].sort()
        if [x[0] for x in v['negatives']]!=list(range(1,len(v['negatives'])+1)): raise ValueError('query rank closure failure')
        v['negatives']=[x[1] for x in v['negatives']]
    return q

def query_weights(queries,qids):
    tree=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for qid in qids:
        z=queries[qid]; tree[z['cluster']][z['structure']][z['anchor']].append(qid)
    ng=len(tree); omega={}
    for g,ss in tree.items():
        for s,aa in ss.items():
            for a,qq in aa.items():
                for qid in qq: omega[qid]=1.0/(ng*len(ss)*len(aa)*len(qq))
    return omega

def training_matrix(queries,qids,emb,return_metadata=False):
    qids=list(qids); omega=query_weights(queries,qids); X=[]; y=[]; w=[]; row_ids=[]
    for qid in qids:
        q=queries[qid]; k=len(q['negatives']); om=omega[qid]
        X.append(pair_feature(emb[q['anchor']],emb[q['true']])); y.append(1.0); w.append(om/2); row_ids.append(qid+'|POS')
        for rank,n in enumerate(q['negatives'],1):
            X.append(pair_feature(emb[q['anchor']],emb[n])); y.append(-1.0); w.append(om/(2*k)); row_ids.append(f'{qid}|NEG|{rank}|{n}')
    X=np.asarray(X,dtype=np.float64); y=np.asarray(y); w=np.asarray(w); w/=w.mean()
    if not return_metadata: return X,y,w
    clusters=sorted({queries[q]['cluster'] for q in qids})
    meta={'training_query_count':len(qids),'training_row_count':len(row_ids),'training_cluster_count':len(clusters),
          'training_query_order_sha256':_hash_lines(qids),'training_cluster_set_sha256':_hash_lines(clusters),
          'training_row_identity_sha256':_hash_lines(row_ids),'normalized_row_weight_sha256':_hash_f64(w),
          'normalized_row_weight_sum':float(w.sum()),'normalized_row_weight_min':float(w.min()),'normalized_row_weight_max':float(w.max())}
    return X,y,w,meta

def fit_logistic(queries,qids,emb,lam,solver):
    X,y,w,provenance=training_matrix(queries,qids,emb,True)
    mu=np.average(X,axis=0,weights=w); var=np.average((X-mu)**2,axis=0,weights=w); scale=np.sqrt(var); scale[scale==0]=1.0
    Z=(X-mu)/scale; d=Z.shape[1]
    def fg(t):
        beta=t[:d]; b=t[d]; margin=y*(Z@beta+b)
        loss=np.sum(w*np.logaddexp(0.0,-margin))+lam*np.dot(beta,beta)
        a=-w*y/(1.0+np.exp(np.clip(margin,-700,700)))
        return float(loss),np.r_[Z.T@a+2*lam*beta,np.sum(a)]
    res=minimize(fg,np.zeros(d+1),jac=True,method='L-BFGS-B',options={'maxiter':solver['max_iter'],'ftol':solver['ftol'],'gtol':solver['gtol'],'maxls':solver['max_line_search']})
    if not res.success: raise RuntimeError('ridge-logistic convergence failure: '+res.message)
    return {'beta':res.x[:d],'intercept':float(res.x[d]),'mean':mu,'scale':scale,'iterations':int(res.nit),'lambda':float(lam),'training_provenance':provenance}

def decision(model,x): return float(np.dot(model['beta'],(x-model['mean'])/model['scale'])+model['intercept'])

def predict_queries(model,queries,qids,emb,tau):
    out={}
    for qid in qids:
        q=queries[qid]; hp=decision(model,pair_feature(emb[q['anchor']],emb[q['true']]))
        hn=[decision(model,pair_feature(emb[q['anchor']],emb[n])) for n in q['negatives']]
        delta=np.asarray([hp-x for x in hn]); credit=np.where(delta>tau,1.0,np.where(delta < -tau,-1.0,0.0))
        out[qid]={'positive_score':hp,'negative_scores':hn,'credits':credit.tolist(),'query_credit':float(credit.mean())}
    return out

def aggregate(queries,pred,qids=None,return_hierarchy=False):
    qids=list(pred) if qids is None else list(qids); tree=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for qid in qids:
        q=queries[qid]; tree[q['cluster']][q['structure']][q['anchor']].append((qid,pred[qid]['query_credit']))
    anchors={}; structures={}; clusters={}
    for g,ss in tree.items():
        sv=[]
        for s,aa in ss.items():
            av=[]
            for a,qq in aa.items():
                v=float(np.mean([x[1] for x in qq])); anchors[f'{g}|{s}|{a}']=v; av.append(v)
            structures[f'{g}|{s}']=float(np.mean(av)); sv.append(structures[f'{g}|{s}'])
        clusters[g]=float(np.mean(sv))
    f=float(np.mean(list(clusters.values())))
    if return_hierarchy: return f,clusters,{'query':{q:pred[q]['query_credit'] for q in qids},'anchor':anchors,'structure':structures,'cluster':clusters}
    return f,clusters

def select_lambda(queries,train_qids,emb,inner_fold,method):
    scores={}; fits={}
    for li,lam in enumerate(method['lambda_grid']):
        oof={}
        for f in range(method['inner_folds']):
            tr=[q for q in train_qids if inner_fold[queries[q]['cluster']]!=f]; te=[q for q in train_qids if inner_fold[queries[q]['cluster']]==f]
            model=fit_logistic(queries,tr,emb,lam,method['solver']); fits[f'lambda_{li}_fold_{f}']=model
            oof.update(predict_queries(model,queries,te,emb,method['tie_tolerance']))
        scores[str(lam)]=aggregate(queries,oof,train_qids)[0]
    best=max(scores.values()); chosen=max(float(k) for k,v in scores.items() if v==best)
    return chosen,scores,fits

def nested_crossfit(queries,emb,split,method):
    pred={}; selection={}; fits={'inner':{},'outer':{}}
    for outer in range(method['outer_folds']):
        train=[q for q,z in queries.items() if split['outer'][z['cluster']]!=outer]; test=[q for q,z in queries.items() if split['outer'][z['cluster']]==outer]
        lam,scores,inner_fits=select_lambda(queries,train,emb,split['inner'][str(outer)],method)
        model=fit_logistic(queries,train,emb,lam,method['solver']); pred.update(predict_queries(model,queries,test,emb,method['tie_tolerance']))
        selection[str(outer)]={'selected_lambda':lam,'inner_scores':scores}; fits['inner'][str(outer)]=inner_fits; fits['outer'][str(outer)]=model
    if set(pred)!=set(queries): raise RuntimeError('OOF query closure failure')
    r,cluster,hierarchy=aggregate(queries,pred,return_hierarchy=True)
    return r,pred,cluster,selection,fits,hierarchy

def bootstrap_cluster_means(cluster_values,draws,seed):
    x=np.asarray([cluster_values[k] for k in sorted(cluster_values)],dtype=np.float64); rng=np.random.Generator(np.random.PCG64(seed)); vals=np.empty(draws)
    for i in range(0,draws,512):
        j=min(draws,i+512); vals[i:j]=x[rng.integers(0,len(x),size=(j-i,len(x)))].mean(axis=1)
    return [float(v) for v in np.percentile(vals,[2.5,97.5])]
