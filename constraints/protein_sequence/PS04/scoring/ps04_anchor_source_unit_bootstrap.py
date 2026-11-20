"""PS04 anchor-source-unit bootstrap.
Multiplicity is used only when source units form HSF values. Once an HSF is
formed, topology and architecture aggregate present biological children equally.
"""
import numpy as np
from collections import defaultdict

def _mean_by_code(values, weights, code, n):
    out=np.zeros((values.shape[0],n),dtype=np.float64); den=np.zeros((values.shape[0],n),dtype=np.float64)
    for g in range(n):
        ix=np.flatnonzero(code==g)
        if len(ix):
            w=weights[:,ix]; den[:,g]=w.sum(axis=1); out[:,g]=np.divide((values[:,ix]*w).sum(axis=1),den[:,g],out=np.zeros(values.shape[0]),where=den[:,g]>0)
    return out,den

def build_maps(anchor_meta, unit_labels):
    # source parent and source unit are immutable point-hierarchy levels.
    def codes(keys):
        vals=sorted(set(keys)); d={x:i for i,x in enumerate(vals)}; return np.array([d[x] for x in keys],dtype=np.int32),vals
    sp,sp_labels=codes([tuple(m[k] for k in ('anchor_hsf_id','anchor_topology_id','anchor_architecture_id','anchor_parent_id')) for m in anchor_meta])
    su,su_labels=codes([tuple(m[k] for k in ('anchor_hsf_id','anchor_topology_id','anchor_architecture_id','anchor_source_unit_id')) for m in anchor_meta])
    # Every source-parent node must lie in one source unit for this release.
    for g in range(len(sp_labels)):
        if len(set(unit_labels[i] for i in np.flatnonzero(sp==g)))!=1: raise ValueError('source parent maps to multiple bootstrap units')
    hsf,hsf_labels=codes([tuple(m[k] for k in ('anchor_architecture_id','anchor_topology_id','anchor_hsf_id')) for m in anchor_meta])
    topo,topo_labels=codes([tuple(m[k] for k in ('anchor_architecture_id','anchor_topology_id')) for m in anchor_meta])
    arch,arch_labels=codes([m['anchor_architecture_id'] for m in anchor_meta])
    # Each source-unit node belongs to one HSF; each HSF to one topology; each topology to one architecture.
    def parent_map(child,parent,nchild):
        out=np.empty(nchild,dtype=np.int32)
        for g in range(nchild):
            q=set(parent[i] for i in np.flatnonzero(child==g))
            if len(q)!=1: raise ValueError('non-nested PS04 hierarchy')
            out[g]=q.pop()
        return out
    # The frozen uncertainty ledger may deliberately merge several point-level
    # source units into one resampling block. It is the authority for sampling.
    unit_code, unit_values = codes(list(unit_labels))
    su_to_bootstrap = parent_map(su, unit_code, len(su_labels))
    return {'anchor_to_sp':sp,'sp_to_su':parent_map(sp,su,len(sp_labels)),'su_to_hsf':parent_map(su,hsf,len(su_labels)),'su_to_bootstrap':su_to_bootstrap,'hsf_to_topo':parent_map(hsf,topo,len(hsf_labels)),'topo_to_arch':parent_map(topo,arch,len(topo_labels)),'n_sp':len(sp_labels),'n_su':len(su_labels),'n_bootstrap':len(unit_values),'n_hsf':len(hsf_labels),'n_topo':len(topo_labels),'n_arch':len(arch_labels)}

def point_source_unit_values(anchor_values,maps):
    v=np.asarray(anchor_values,dtype=np.float64)[None,:]; one=np.ones_like(v)
    sp,_=_mean_by_code(v,one,maps['anchor_to_sp'],maps['n_sp'])
    su,_=_mean_by_code(sp,np.ones_like(sp),maps['sp_to_su'],maps['n_su'])
    return su[0]

def vectorized_draws(source_unit_values,maps,counts):
    v=np.broadcast_to(np.asarray(source_unit_values,dtype=np.float64),(counts.shape[0],len(source_unit_values)))
    # Unit multiplicity determines the HSF average.
    if counts.shape[1] != maps['n_bootstrap']: raise ValueError('frozen bootstrap count width mismatch')
    source_unit_weights=counts[:,maps['su_to_bootstrap']]
    hsf,den_hsf=_mean_by_code(v,source_unit_weights,maps['su_to_hsf'],maps['n_hsf'])
    # HSF presence, not source-unit multiplicity, is the topology child weight.
    topo,den_topo=_mean_by_code(hsf,(den_hsf>0).astype(float),maps['hsf_to_topo'],maps['n_topo'])
    arch,den_arch=_mean_by_code(topo,(den_topo>0).astype(float),maps['topo_to_arch'],maps['n_arch'])
    overall,den_overall=_mean_by_code(arch,(den_arch>0).astype(float),np.zeros(maps['n_arch'],dtype=np.int32),1)
    return overall[:,0]

def bootstrap(source_unit_values,maps,draws,seed,batch=256):
    n=maps['n_bootstrap']; rng=np.random.Generator(np.random.PCG64(seed)); out=[]
    for start in range(0,draws,batch):
        b=min(batch,draws-start); counts=rng.multinomial(n,np.full(n,1/n),size=b).astype(float)
        out.append(vectorized_draws(source_unit_values,maps,counts))
    return np.concatenate(out)

def reference_draw(source_unit_values,maps,count_row):
    # Independent literal copy construction used only for parity fixtures.
    hsf=defaultdict(list)
    for source_unit, value in enumerate(source_unit_values):
        c=count_row[maps['su_to_bootstrap'][source_unit]]
        for _ in range(int(c)): hsf[maps['su_to_hsf'][source_unit]].append(float(value))
    topo=defaultdict(list)
    for h,vs in hsf.items(): topo[maps['hsf_to_topo'][h]].append(sum(vs)/len(vs))
    arch=defaultdict(list)
    for t,vs in topo.items(): arch[maps['topo_to_arch'][t]].append(sum(vs)/len(vs))
    vals=[sum(vs)/len(vs) for vs in arch.values()]
    return sum(vals)/len(vals)

def parity_fixture():
    # U1 copied twice must not upweight its topology after HSF formation.
    maps={'su_to_hsf':np.array([0,1]),'su_to_bootstrap':np.array([0,1]),'hsf_to_topo':np.array([0,1]),'topo_to_arch':np.array([0,0]),'n_bootstrap':2,'n_hsf':2,'n_topo':2,'n_arch':1}
    values=np.array([0.,1.]); counts=np.array([[2.,1.],[3.,1.],[1.,2.]])
    got=vectorized_draws(values,maps,counts); want=np.array([reference_draw(values,maps,x) for x in counts])
    if not np.array_equal(got,want): raise AssertionError((got,want))
    if not np.array_equal(got,np.array([.5,.5,.5])): raise AssertionError(got)
    return {'draws':3,'max_abs_difference':float(np.max(np.abs(got-want))),'values':got.tolist()}
