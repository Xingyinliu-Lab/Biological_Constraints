from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class WindowProjection:
    start:int; end:int; positions:np.ndarray; token_coverage:np.ndarray

def frozen_windows(length):
    length=int(length)
    if length<=0:raise ValueError('empty sequence')
    if length<=1000:return ((0,length),)
    starts=list(range(0,length-1000+1,500));last=length-1000
    if starts[-1]!=last:starts.append(last)
    return tuple((s,s+1000) for s in starts)

def native_single_base_spans(length):return tuple((i,i+1) for i in range(length))

def project_tokens_to_nucleotides(hidden,spans,nucleotide_length):
    h=np.asarray(hidden,dtype=np.float64)
    if h.ndim!=2 or h.shape[0]!=len(spans) or not np.isfinite(h).all():raise ValueError('invalid token states')
    sums=np.zeros((nucleotide_length,h.shape[1]));counts=np.zeros(nucleotide_length,dtype=np.int64)
    for row,(a,b) in zip(h,spans):
        a,b=int(a),int(b)
        if a==b==0:continue
        if not 0<=a<b<=nucleotide_length:raise ValueError('invalid token span')
        sums[a:b]+=row;counts[a:b]+=1
    if np.any(counts==0):raise ValueError('incomplete nucleotide coverage')
    return (sums/counts[:,None]).astype(np.float32),counts

def merge_window_projections(length,projections):
    projections=list(projections);expected=frozen_windows(length);observed=tuple((int(x.start),int(x.end)) for x in projections)
    if observed!=expected:raise ValueError(f'window mismatch: {observed} != {expected}')
    sums=None;counts=np.zeros(length,dtype=np.int64)
    for p in projections:
        arr=np.asarray(p.positions,dtype=np.float64)
        if arr.ndim!=2 or arr.shape[0]!=p.end-p.start or not np.isfinite(arr).all() or np.any(np.asarray(p.token_coverage)<=0):raise ValueError('invalid window projection')
        if sums is None:sums=np.zeros((length,arr.shape[1]))
        sums[p.start:p.end]+=arr;counts[p.start:p.end]+=1
    if sums is None or np.any(counts==0):raise ValueError('incomplete global coverage')
    return (sums/counts[:,None]).astype(np.float32),counts

def pool_codon_local(position_states,core_start,core_end):
    x=np.asarray(position_states,dtype=np.float64);core_start=int(core_start);core_end=int(core_end)
    if x.ndim!=2 or not np.isfinite(x).all() or core_end-core_start!=3 or not 0<=core_start<core_end<=len(x):raise ValueError('invalid three-nucleotide codon interval')
    return x[core_start:core_end].mean(0).astype(np.float32)
