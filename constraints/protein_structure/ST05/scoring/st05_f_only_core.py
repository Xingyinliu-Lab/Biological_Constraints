from __future__ import annotations

from collections import defaultdict

import numpy as np


BOOTSTRAP_DRAWS=20000
BOOTSTRAP_SEED=15005


def _mean(values: list[float]) -> float:
    if not values:
        raise ValueError("empty aggregation unit")
    return float(np.mean(np.asarray(values,dtype=np.float64)))


def aggregate(rows: list[dict[str,object]]) -> dict[str,object]:
    by_base: dict[str,list[float]]=defaultdict(list)
    base_meta: dict[str,tuple[str,str]]={}
    for row in rows:
        base=str(row["base_object_id"]); parent=str(row["source_parent_id"]); cluster=str(row["source_dependency_cluster_id"])
        credit=float(row["F_credit"])
        if not np.isfinite(credit) or credit not in {0.0,1.0}:
            raise ValueError("invalid binary satisfaction credit")
        if base in base_meta and base_meta[base]!=(parent,cluster):
            raise ValueError("base hierarchy conflict")
        base_meta[base]=(parent,cluster); by_base[base].append(credit)
    base_scores={base:_mean(values) for base,values in by_base.items()}
    by_parent: dict[tuple[str,str],list[float]]=defaultdict(list)
    for base,score in base_scores.items():
        parent,cluster=base_meta[base]; by_parent[(cluster,parent)].append(score)
    parent_scores={key:_mean(values) for key,values in by_parent.items()}
    by_cluster: dict[str,list[float]]=defaultdict(list)
    for (cluster,_parent),score in parent_scores.items(): by_cluster[cluster].append(score)
    cluster_scores={cluster:_mean(values) for cluster,values in by_cluster.items()}
    return {
        "F":_mean(list(cluster_scores.values())),
        "base_scores":base_scores,
        "parent_scores":{f"{cluster}\t{parent}":score for (cluster,parent),score in parent_scores.items()},
        "cluster_scores":cluster_scores,
    }


def bootstrap_interval(cluster_scores: dict[str,float],draws: int=BOOTSTRAP_DRAWS,seed: int=BOOTSTRAP_SEED) -> tuple[float,float]:
    cluster_ids=sorted(cluster_scores)
    values=np.asarray([cluster_scores[x] for x in cluster_ids],dtype=np.float64)
    if len(values)==0: raise ValueError("no dependency clusters")
    rng=np.random.default_rng(seed)
    estimates=np.empty(draws,dtype=np.float64)
    for draw in range(draws):
        indices=rng.integers(0,len(values),size=len(values))
        estimates[draw]=float(np.mean(values[indices]))
    lo,hi=np.percentile(estimates,[2.5,97.5])
    return float(lo),float(hi)
