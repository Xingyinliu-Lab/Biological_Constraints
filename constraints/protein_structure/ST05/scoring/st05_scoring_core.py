from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_core(path: Path):
    spec=importlib.util.spec_from_file_location("st05_frozen_numeric_core",path)
    if spec is None or spec.loader is None: raise ValueError("numeric core import failure")
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if Path(module.__file__).resolve()!=path.resolve(): raise ValueError("numeric core path mismatch")
    return module


def string_array(values) -> list[str]:
    return [str(x) for x in np.asarray(values).tolist()]


def weighted_norm(x: np.ndarray,weights: np.ndarray|None) -> float:
    flat=np.asarray(x,dtype=np.float64).reshape(-1)
    if not np.isfinite(flat).all(): raise ValueError("nonfinite representation/residual")
    maximum=float(np.max(np.abs(flat))) if flat.size else 0.0
    if maximum==0.0: return 0.0
    scaled=flat/maximum
    if weights is None:
        value=maximum*float(np.sqrt(np.sum(scaled*scaled,dtype=np.float64)))
        if not np.isfinite(value): raise ValueError("nonfinite norm result")
        return value
    if weights.shape!=flat.shape or not np.isfinite(weights).all() or np.any(weights<=0):
        raise ValueError("invalid frozen component weights")
    value=maximum*float(np.sqrt(np.sum(weights*scaled*scaled,dtype=np.float64)))
    if not np.isfinite(value): raise ValueError("nonfinite weighted norm result")
    return value


def evaluate_law(z: np.ndarray,zg: np.ndarray,matrix: np.ndarray,bias: np.ndarray,weights: np.ndarray|None,tau_abs: float,tau_rel: float) -> dict[str,float]:
    expected=matrix@z+bias
    if not np.isfinite(expected).all(): raise ValueError("nonfinite declared-action output")
    residual=weighted_norm(zg-expected,weights)
    scale=max(weighted_norm(zg,weights),weighted_norm(expected,weights))
    tau=tau_abs+tau_rel*scale
    if not np.isfinite([residual,scale,tau]).all(): raise ValueError("nonfinite residual/scale/tolerance")
    ratio=residual/(tau+np.finfo(np.float64).eps)
    if not np.isfinite(ratio): raise ValueError("nonfinite residual/tolerance ratio")
    return {"residual":residual,"scale":scale,"tau":tau,"ratio":ratio,"credit":float(residual<=tau)}


def validate_operator_evidence(evidence_path: Path,report: dict,dimension: int,tolerance: float) -> None:
    if report.get("status")!="PASS": raise ValueError("operator gate status failure")
    if report.get("operator_evidence_sha256")!=sha256(evidence_path): raise ValueError("operator evidence hash mismatch")
    evidence=np.load(evidence_path,allow_pickle=False)
    identity_matrix=np.asarray(evidence["identity_matrix"],dtype=np.float64); identity_bias=np.asarray(evidence["identity_bias"],dtype=np.float64)
    if identity_matrix.shape!=(dimension,dimension) or identity_bias.shape!=(dimension,): raise ValueError("operator identity evidence shape failure")
    errors=[float(np.max(np.abs(identity_matrix-np.eye(dimension)))),float(np.max(np.abs(identity_bias)))]
    ca=np.asarray(evidence["composition_a_matrices"],dtype=np.float64); cb=np.asarray(evidence["composition_b_matrices"],dtype=np.float64); cc=np.asarray(evidence["composition_declared_matrices"],dtype=np.float64)
    ba=np.asarray(evidence["composition_a_biases"],dtype=np.float64); bb=np.asarray(evidence["composition_b_biases"],dtype=np.float64); bc=np.asarray(evidence["composition_declared_biases"],dtype=np.float64)
    if ca.ndim!=3 or ca.shape!=cb.shape or ca.shape!=cc.shape or ca.shape[1:]!=(dimension,dimension) or ba.shape!=bb.shape or ba.shape!=bc.shape or ba.shape!=(ca.shape[0],dimension) or ca.shape[0]==0: raise ValueError("operator composition evidence shape failure")
    for i in range(ca.shape[0]):
        errors.extend([float(np.max(np.abs(ca[i]@cb[i]-cc[i]))),float(np.max(np.abs(ca[i]@bb[i]+ba[i]-bc[i])))])
    fm=np.asarray(evidence["inverse_forward_matrices"],dtype=np.float64); fb=np.asarray(evidence["inverse_forward_biases"],dtype=np.float64); im=np.asarray(evidence["inverse_declared_matrices"],dtype=np.float64); ib=np.asarray(evidence["inverse_declared_biases"],dtype=np.float64)
    if fm.ndim!=3 or fm.shape!=im.shape or fm.shape[1:]!=(dimension,dimension) or fb.shape!=ib.shape or fb.shape!=(fm.shape[0],dimension) or fm.shape[0]==0: raise ValueError("operator inverse evidence shape failure")
    for i in range(fm.shape[0]):
        errors.extend([float(np.max(np.abs(im[i]@fm[i]-np.eye(dimension)))),float(np.max(np.abs(im[i]@fb[i]+ib[i])))])
    if not np.isfinite(errors).all() or max(errors)>tolerance: raise ValueError("operator behavior gate failure")
