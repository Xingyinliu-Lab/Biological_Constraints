from __future__ import annotations
from collections import defaultdict
import numpy as np

def cosine_distance(a, b):
    a=np.asarray(a,dtype=np.float64);b=np.asarray(b,dtype=np.float64)
    na=float(np.linalg.norm(a));nb=float(np.linalg.norm(b))
    if not np.isfinite(a).all() or not np.isfinite(b).all() or na<=0 or nb<=0: raise ValueError("invalid embedding")
    return float(1.0-np.dot(a,b)/(na*nb))

def tie_F_credit(d_positive, d_negative, tau):
    if not all(np.isfinite(x) for x in (d_positive,d_negative,tau)) or tau<0: raise ValueError("invalid distance/tau")
    delta=d_negative-d_positive
    return 1.0 if delta>tau else -1.0 if delta < -tau else 0.0

def aggregate_query_rows(rows):
    by_cluster_transition=defaultdict(list)
    for r in rows: by_cluster_transition[(r["dependency_cluster_id"],r["transition_id"])].append(float(r["query_F"]))
    cluster_transition={k:float(np.mean(v)) for k,v in by_cluster_transition.items()}
    by_transition=defaultdict(list)
    for (cluster,transition),score in cluster_transition.items(): by_transition[transition].append(score)
    transition={k:float(np.mean(v)) for k,v in by_transition.items()}
    return float(np.mean(list(transition.values()))),transition,cluster_transition

def cluster_bootstrap_distribution(rows, draws=20000, seed=20261014):
    clusters=sorted({r["dependency_cluster_id"] for r in rows})
    transitions=sorted({r["transition_id"] for r in rows})
    ci={value:index for index,value in enumerate(clusters)};ti={value:index for index,value in enumerate(transitions)}
    grouped=defaultdict(list)
    for r in rows:grouped[(r["dependency_cluster_id"],r["transition_id"])].append(float(r["query_F"]))
    matrix=np.full((len(clusters),len(transitions)),np.nan,dtype=np.float64)
    for (cluster,transition),values in grouped.items():matrix[ci[cluster],ti[transition]]=float(np.mean(values))
    rng=np.random.default_rng(seed);values=np.empty(draws,dtype=np.float64)
    for draw in range(draws):
        sample=rng.choice(len(clusters),size=len(clusters),replace=True)
        values[draw]=float(np.mean(np.nanmean(matrix[sample],axis=0)))
    return values


def cluster_bootstrap(rows, draws=20000, seed=20261014):
    values = cluster_bootstrap_distribution(rows, draws, seed)
    return [float(np.quantile(values,.025)),float(np.quantile(values,.975))]
