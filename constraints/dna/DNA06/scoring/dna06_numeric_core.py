from __future__ import annotations
import hashlib, math
import numpy as np

LAMBDA_GRID=(1e-4,1e-3,1e-2,1e-1,1.0,10.0,100.0,1000.0,10000.0)

def fold_map(cluster_ids, folds, seed):
    keys=sorted(set(cluster_ids), key=lambda x: hashlib.sha256((seed+'|'+x).encode()).hexdigest())
    return {x:i%folds for i,x in enumerate(keys)}

def standardize_fit(x):
    mu=x.mean(0); sd=x.std(0); sd=np.where(sd>1e-12,sd,1.0)
    return mu,sd

def ridge_logistic_fit(x,y,lam,max_iter=1000,tol=1e-8):
    x=np.asarray(x,float); y=np.asarray(y,float)
    if x.ndim!=2 or len(x)!=len(y) or not np.isfinite(x).all() or not np.isfinite(y).all() or not math.isfinite(lam) or lam<=0 or set(np.unique(y))!={0.0,1.0}: raise ValueError('invalid binary training data')
    from scipy.optimize import minimize
    xx=np.column_stack([np.ones(len(x)),x]); pen=np.r_[0.0,np.full(x.shape[1],lam)]
    def objective(beta):
        eta=xx@beta
        value=np.logaddexp(0.0,eta).sum()-y@eta+0.5*np.dot(pen*beta,beta)
        prob=np.empty_like(eta); pos=eta>=0
        prob[pos]=1.0/(1.0+np.exp(-eta[pos])); exp_eta=np.exp(eta[~pos]); prob[~pos]=exp_eta/(1.0+exp_eta)
        grad=xx.T@(prob-y)+pen*beta
        if not math.isfinite(float(value)) or not np.isfinite(grad).all(): raise FloatingPointError('non-finite ridge objective')
        return float(value),grad
    fit=minimize(objective,np.zeros(xx.shape[1]),jac=True,method='L-BFGS-B',options={'maxiter':max_iter,'ftol':tol,'gtol':tol,'maxls':50})
    if not fit.success: raise RuntimeError('ridge solver did not converge: '+str(fit.message))
    beta=np.asarray(fit.x,float)
    if not np.isfinite(beta).all(): raise FloatingPointError('non-finite ridge solution')
    return beta

def predict(x,beta):
    x=np.asarray(x,float); beta=np.asarray(beta,float)
    if not np.isfinite(x).all() or not np.isfinite(beta).all(): raise ValueError('non-finite prediction input')
    out=np.column_stack([np.ones(len(x)),x])@beta
    if not np.isfinite(out).all(): raise FloatingPointError('non-finite prediction output')
    return out

def pair_F_credit(m,s,tau=1e-8):
    m,s,tau=float(m),float(s),float(tau)
    if not all(map(math.isfinite,(m,s,tau))) or tau<0: raise ValueError('non-finite score or invalid tie tolerance')
    d=m-s
    return 1.0 if d>tau else (-1.0 if d<-tau else 0.0)

def site_scores(rows,tau=1e-8):
    by={}
    for r in rows: by.setdefault(r['site_id'],[]).append(r)
    out=[]
    for site,rr in by.items():
        if len({(x['parent_cds_id'],x['dependency_cluster_id']) for x in rr})!=1: raise ValueError('site identity conflict: '+site)
        labels={x['label'] for x in rr}
        if not labels <= {'SYNONYMOUS','ORDINARY_MISSENSE'}: raise ValueError('invalid consequence label')
        mis=[float(x['score']) for x in rr if x['label']=='ORDINARY_MISSENSE']; syn=[float(x['score']) for x in rr if x['label']=='SYNONYMOUS']
        if not mis or not syn: raise ValueError('site missing class: '+site)
        val=sum(pair_F_credit(m,s,tau) for m in mis for s in syn)/(len(mis)*len(syn))
        if not math.isfinite(val) or not -1<=val<=1: raise ValueError('invalid site score')
        q=rr[0]; out.append({'site_id':site,'parent_cds_id':q['parent_cds_id'],'dependency_cluster_id':q['dependency_cluster_id'],'site_F':val})
    return out

def pair_credit_ledger(rows,tau=1e-8):
    """Materialize frozen missense-versus-synonymous credits per site."""
    by={}
    for row in rows: by.setdefault(row['site_id'],[]).append(row)
    out=[]
    for site_id,site_rows in by.items():
        identity={(row['parent_cds_id'],row['dependency_cluster_id']) for row in site_rows}
        if len(identity)!=1: raise ValueError('site identity conflict: '+site_id)
        missense=[row for row in site_rows if row['label']=='ORDINARY_MISSENSE']; synonymous=[row for row in site_rows if row['label']=='SYNONYMOUS']
        if not missense or not synonymous or len(missense)+len(synonymous)!=len(site_rows): raise ValueError('site missing or invalid class: '+site_id)
        for missense_row in missense:
            for synonymous_row in synonymous:
                missense_score=float(missense_row['score']); synonymous_score=float(synonymous_row['score']); margin=missense_score-synonymous_score
                out.append({'pair_id':f"{site_id}|{missense_row['variant_id']}|{synonymous_row['variant_id']}",'site_id':site_id,'parent_cds_id':missense_row['parent_cds_id'],'dependency_cluster_id':missense_row['dependency_cluster_id'],'missense_variant_id':missense_row['variant_id'],'synonymous_variant_id':synonymous_row['variant_id'],'missense_decision_score':missense_score,'synonymous_decision_score':synonymous_score,'margin':margin,'tie':abs(margin)<=tau,'F_credit':pair_F_credit(missense_score,synonymous_score,tau)})
    if len({row['pair_id'] for row in out})!=len(out): raise ValueError('duplicate matched-pair identity')
    return out

def site_scores_from_pair_ledger(pair_rows):
    grouped={}
    for row in pair_rows: grouped.setdefault(row['site_id'],[]).append(row)
    out=[]
    for site_id,rows in grouped.items():
        identity={(row['parent_cds_id'],row['dependency_cluster_id']) for row in rows}
        if len(identity)!=1: raise ValueError('pair ledger site identity conflict: '+site_id)
        parent_id,cluster_id=next(iter(identity)); out.append({'site_id':site_id,'parent_cds_id':parent_id,'dependency_cluster_id':cluster_id,'site_F':sum(float(row['F_credit']) for row in rows)/len(rows)})
    return out

def parent_scores(site_rows):
    grouped={}
    for row in site_rows: grouped.setdefault((row['dependency_cluster_id'],row['parent_cds_id']),[]).append(float(row['site_F']))
    return [{'parent_cds_id':parent_id,'dependency_cluster_id':cluster_id,'site_count':len(values),'parent_F':sum(values)/len(values)} for (cluster_id,parent_id),values in sorted(grouped.items())]

def hierarchical(site_rows):
    clusters={}
    for row in parent_scores(site_rows): clusters.setdefault(row['dependency_cluster_id'],[]).append(float(row['parent_F']))
    cr={c:sum(v)/len(v) for c,v in clusters.items()}
    if not cr or not np.isfinite(list(cr.values())).all() or any(not -1<=v<=1 for v in cr.values()): raise ValueError('invalid cluster scores')
    return sum(cr.values())/len(cr),cr

def bootstrap_distribution(cluster_scores,draws=20000,seed=20261014):
    ids=sorted(cluster_scores); vals=np.array([cluster_scores[x] for x in ids],float); rng=np.random.default_rng(seed)
    sims=np.empty(draws)
    for i in range(draws): sims[i]=vals[rng.integers(0,len(vals),len(vals))].mean()
    if not np.isfinite(sims).all(): raise FloatingPointError('non-finite bootstrap')
    return sims

def bootstrap(cluster_scores,draws=20000,seed=20261014):
    sims=bootstrap_distribution(cluster_scores,draws,seed)
    return [float(np.quantile(sims,.025)),float(np.quantile(sims,.975))]

def inner_select(x,y,clusters,sites,parents):
    fm=fold_map(clusters,4,'DNA06_INNER_V010'); perf={}
    for lam in LAMBDA_GRID:
        pred=np.empty(len(y))
        for f in range(4):
            tr=np.array([fm[c]!=f for c in clusters]); va=~tr; mu,sd=standardize_fit(x[tr]); b=ridge_logistic_fit((x[tr]-mu)/sd,y[tr],lam); pred[va]=predict((x[va]-mu)/sd,b)
        rr=[{'site_id':sites[i],'parent_cds_id':parents[i],'dependency_cluster_id':clusters[i],'label':'ORDINARY_MISSENSE' if y[i]==1 else 'SYNONYMOUS','score':pred[i]} for i in range(len(y))]
        perf[lam]=hierarchical(site_scores(rr))[0]
    best=max(perf.values()); return max(k for k,v in perf.items() if abs(v-best)<=1e-12),perf
