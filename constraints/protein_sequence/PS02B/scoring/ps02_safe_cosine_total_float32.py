"""Canonical PS02 total float32 safe-cosine operator."""
import math
import numpy as np
OPERATOR_ID='PS02_SAFE_COSINE_TOTAL_FLOAT32'
def safe_cosine_total_float32(a,b):
    a=np.asarray(a,dtype=np.float32); b=np.asarray(b,dtype=np.float32)
    if a.ndim!=1 or b.ndim!=1 or a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('invalid finite equal-length representation vectors')
    na=np.linalg.norm(a); nb=np.linalg.norm(b)
    if na==0 or nb==0:return 0.0
    v=float(np.dot(a,b)/(na*nb))
    if not math.isfinite(v):raise ValueError('nonfinite safe cosine')
    return v