"""PS05 formal uncertainty kernel: exact hierarchy, vectorized nested bootstrap."""
import importlib.util
from collections import defaultdict
import numpy as np
SEED=20261008
DRAWS=20000
DELTA_NUM=0.0
def load_scorer(path):
 s=importlib.util.spec_from_file_location('ps05_formal',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def bootstrap(rows,population_identity,scorer_path):
 m=load_scorer(scorer_path); trace,_=m.score(rows,population_identity)
 groups=defaultdict(lambda:defaultdict(list)); topo_arch={}
 for n in trace['SOURCE_PARENT']:
  meta=n['meta'];topo=meta['anchor_topology_id'];hsf=meta['anchor_hsf_id'];arch=meta['anchor_architecture_id'];groups[topo][hsf].append(float(n['value']));topo_arch[topo]=arch
 rng=np.random.Generator(np.random.PCG64(SEED));topo_samples={}
 for topo,hsfs in groups.items():
  keys=sorted(hsfs); n_h=len(keys); chosen=rng.integers(n_h,size=(DRAWS,n_h)); vals=np.empty((DRAWS,n_h),dtype=np.float64)
  for j in range(n_h):
   chosen_h=chosen[:,j]; col=np.empty(DRAWS,dtype=np.float64)
   for hi,h in enumerate(keys):
    mask=chosen_h==hi; count=int(mask.sum()); parents=np.asarray(hsfs[h],dtype=np.float64); n_p=len(parents)
    col[mask]=parents[rng.integers(n_p,size=(count,n_p))].mean(axis=1)
   vals[:,j]=col
  topo_samples[topo]=vals.mean(axis=1)
 arches=defaultdict(list)
 for topo,v in topo_samples.items():arches[topo_arch[topo]].append(v)
 arch_values=[np.mean(vs,axis=0) for _,vs in sorted(arches.items())]
 return np.mean(arch_values,axis=0)
def ci95(values):
 return tuple(float(x) for x in np.quantile(values,[0.025,0.975],method='linear'))
