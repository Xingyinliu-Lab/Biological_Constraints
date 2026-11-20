"""Shared total float32 cosine operator for Biological Constraints."""
import numpy as np
OPERATOR_ID = "BC_SAFE_COSINE_TOTAL_FLOAT32_v1"
def safe_cosine_total_float32(a,b):
    x=np.asarray(a,dtype=np.float32); y=np.asarray(b,dtype=np.float32)
    if x.ndim!=1 or y.ndim!=1 or x.shape!=y.shape: raise ValueError("equal-length vectors required")
    if not np.isfinite(x).all() or not np.isfinite(y).all(): raise ValueError("non-finite vector")
    nx=np.linalg.norm(x); ny=np.linalg.norm(y)
    return np.float32(0.0) if nx==0 or ny==0 else np.float32(np.dot(x,y)/(nx*ny))