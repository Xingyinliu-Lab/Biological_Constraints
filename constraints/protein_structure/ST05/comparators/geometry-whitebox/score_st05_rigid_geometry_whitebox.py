from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scoring"))
from st05_f_only_core import aggregate


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()


def read_pdb(path: Path) -> tuple[np.ndarray,str]:
    coords=[]; identity=[]
    for line in path.read_text(encoding="ascii",errors="strict").splitlines(keepends=True):
        if line[:6] not in {"ATOM  ","HETATM"}: continue
        if len(line)<54: raise ValueError(f"short coordinate record: {path}")
        xyz=np.asarray([float(line[30:38]),float(line[38:46]),float(line[46:54])],dtype=np.float64)
        if not np.isfinite(xyz).all(): raise ValueError(f"nonfinite coordinates: {path}")
        coords.append(xyz); identity.append(line[:30]+line[54:])
    if not coords: raise ValueError(f"no coordinates: {path}")
    return np.asarray(coords,dtype=np.float64),hashlib.sha256("".join(identity).encode("ascii")).hexdigest()


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--release",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    args=ap.parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True)
    tau_abs=float(np.sqrt(3*4680)*0.000501)
    tau_rel=float(64*np.finfo(np.float64).eps)
    with (args.release/"objects.tsv").open(encoding="utf-8",newline="") as h: objects={r["base_object_id"]:r for r in csv.DictReader(h,delimiter="\t")}
    with (args.release/"relations.tsv").open(encoding="utf-8",newline="") as h: relations=list(csv.DictReader(h,delimiter="\t"))
    if len(objects)!=60 or len(relations)!=960: raise ValueError("population count failure")
    if not np.isfinite([tau_abs,tau_rel]).all() or tau_abs<0 or tau_rel<0: raise ValueError("invalid preflight tolerance")
    rows=[]; by_base=defaultdict(set); residuals=[]; ratios=[]; scope_checks=0
    for relation in relations:
        obj=objects[relation["base_object_id"]]
        base=args.release/obj["coordinate_locator"]; variant=args.release/relation["variant_coordinate_locator"]
        if sha256(base)!=relation["base_coordinate_sha256"] or sha256(variant)!=relation["variant_coordinate_sha256"]: raise ValueError("coordinate hash failure")
        x,base_scope=read_pdb(base); xg,variant_scope=read_pdb(variant)
        if x.shape!=xg.shape or base_scope!=variant_scope or base_scope!=relation["coordinate_scope_sha256"]: raise ValueError("coordinate scope/order failure")
        matrix=np.asarray(json.loads(relation["rotation_matrix"]),dtype=np.float64); translation=np.asarray(json.loads(relation["translation_vector"]),dtype=np.float64)
        if matrix.shape!=(3,3) or translation.shape!=(3,) or not np.isfinite(matrix).all() or not np.isfinite(translation).all(): raise ValueError("transform format failure")
        if np.max(np.abs(matrix.T@matrix-np.eye(3)))>1e-12 or abs(float(np.linalg.det(matrix))-1.0)>1e-12: raise ValueError("proper rotation failure")
        expected=x@matrix.T+translation
        residual=float(np.linalg.norm(xg-expected)); scale=max(float(np.linalg.norm(xg)),float(np.linalg.norm(expected)))
        tau=tau_abs+tau_rel*scale; ratio=residual/(tau+np.finfo(np.float64).eps); f_credit=float(residual<=tau)
        rows.append({**relation,"residual_frobenius":residual,"scale":scale,"tau":tau,"normalized_residual":ratio,"margin":tau-residual,"F_credit":f_credit,"status":"PASS" if f_credit else "NUMERICAL_REFERENCE_QC_FAILURE","threshold_rule":"residual_frobenius <= tau"})
        residuals.append(residual); ratios.append(ratio); by_base[relation["base_object_id"]].add(relation["transform_id"]); scope_checks+=1
    expected_transforms={r["transform_id"] for r in relations}
    if len(expected_transforms)!=16 or any(value!=expected_transforms for value in by_base.values()): raise ValueError("complete panel failure")
    aggregation=aggregate(rows); score=aggregation["F"]; clusters=aggregation["cluster_scores"]
    status="PASS" if score==1.0 else "NUMERICAL_REFERENCE_QC_FAILURE"
    fields=list(relations[0])+["residual_frobenius","scale","tau","normalized_residual","margin","F_credit","status","threshold_rule"]
    with (args.out_dir/"WB_TRANSFORM_LEDGER.tsv").open("w",encoding="utf-8",newline="") as h:
        writer=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    result={"artifact_id":"ST05_RIGID_GEOMETRY_WHITEBOX_F_ONLY_RESULT_v0.1.2","status":status,"role":"LEVEL_0_ANALYTICAL_WHITEBOX_REFERENCE","formal_interface":"F_ONLY","F":score,"F_95_CI":None,"counts":{"bases":len(objects),"relations":len(relations),"dependency_clusters":len(clusters)},"hierarchy_F":{"base_F":aggregation["base_scores"],"parent_F":aggregation["parent_scores"],"dependency_cluster_F":clusters},"residual_summary":{"min":float(np.min(residuals)),"median":float(np.median(residuals)),"max":float(np.max(residuals))},"normalized_residual_summary":{"min":float(np.min(ratios)),"median":float(np.median(ratios)),"max":float(np.max(ratios))},"current_scorer_sha256":sha256(Path(__file__).resolve()),"ledger_sha256":sha256(args.out_dir/"WB_TRANSFORM_LEDGER.tsv")}
    (args.out_dir/"WB_RESULT.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    validation={"artifact_id":"ST05_RIGID_GEOMETRY_WHITEBOX_F_ONLY_VALIDATION_v0.1.2","status":status,"coordinate_hash_checks":2*len(relations),"coordinate_scope_order_checks":scope_checks,"proper_rotation_checks":len(relations),"complete_panel_checks":len(objects),"nonfinite_values":0,"failed_relations":sum(row["F_credit"]==0 for row in rows)}
    (args.out_dir/"WB_VALIDATION.json").write_text(json.dumps(validation,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__": main()


