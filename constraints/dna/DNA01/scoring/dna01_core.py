#!/usr/bin/env python3
"""Shared frozen utilities for the DNA01 public scoring package."""

import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np


FORMAL_ALPHA = 0.001
FORMAL_DRAWS = 20000
FORMAL_SEED = 20261013
FORMAL_COVERAGE_THRESHOLD = 0.8
FORMAL_TARGET_SCALE_FLOOR = 0.0


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_vectors(path):
    with np.load(path, allow_pickle=False) as payload:
        if set(payload.files) != {"object_ids", "vectors"}:
            raise ValueError("embedding NPZ must contain exactly object_ids and vectors")
        object_ids = payload["object_ids"].astype(str).tolist()
        vectors = np.asarray(payload["vectors"])
    if len(object_ids) != len(set(object_ids)):
        raise ValueError("duplicate embedding object_id")
    if vectors.ndim != 2 or vectors.shape[0] != len(object_ids) or vectors.shape[1] <= 0:
        raise ValueError("embedding vectors must have shape (n, d) with d > 0")
    if not (np.issubdtype(vectors.dtype, np.floating) or np.issubdtype(vectors.dtype, np.integer)):
        raise ValueError("embedding vectors must use a real floating-point or integer dtype")
    if not np.isfinite(vectors).all():
        raise ValueError("embedding vectors must be finite")
    return dict(zip(object_ids, vectors)), vectors.dtype.str


def validate_refset(objects, relations, anchors, folds, dependency_rows, repeat_rows):
    object_by_id = {}
    for row in objects:
        oid = row["object_id"]
        if oid in object_by_id:
            raise ValueError("duplicate object_id")
        object_by_id[oid] = row

    fold_by_cluster = {}
    for row in folds:
        cluster = row["source_dependency_cluster_id"]
        fold = int(row["outer_fold_id"])
        if cluster in fold_by_cluster or fold not in range(5):
            raise ValueError("invalid or duplicate fold assignment")
        fold_by_cluster[cluster] = fold
    if set(fold_by_cluster.values()) != set(range(5)):
        raise ValueError("formal fold assignment must contain folds 0 through 4")

    relation_by_id = {}
    for row in relations:
        rid = row["relation_id"]
        if rid in relation_by_id:
            raise ValueError("duplicate relation_id")
        fold = int(row["outer_fold_id"])
        cluster = row["source_dependency_cluster_id"]
        if fold_by_cluster.get(cluster) != fold:
            raise ValueError("relation fold mismatch")
        fwd = object_by_id.get(row["forward_object_id"])
        rc = object_by_id.get(row["rc_object_id"])
        if fwd is None or rc is None:
            raise ValueError("relation object absent from object table")
        relation_by_id[rid] = row

    anchors_by_relation = defaultdict(list)
    anchor_keys = set()
    for row in anchors:
        key = (row["relation_id"], int(row["anchor_rank"]))
        if key in anchor_keys:
            raise ValueError("duplicate anchor key")
        anchor_keys.add(key)
        if row["relation_id"] not in relation_by_id:
            raise ValueError("anchor references an unknown relation")
        anchors_by_relation[row["relation_id"]].append(row)
    if set(anchors_by_relation) != set(relation_by_id):
        raise ValueError("anchor/relation universe mismatch")

    for rid, selected in anchors_by_relation.items():
        relation = relation_by_id[rid]
        selected.sort(key=lambda row: int(row["anchor_rank"]))
        if len(selected) != 8 or [int(row["anchor_rank"]) for row in selected] != list(range(1, 9)):
            raise ValueError("each relation must have anchor ranks 1 through 8")
        if len({row["anchor_target_object_id"] for row in selected}) != 8:
            raise ValueError("anchor targets must be distinct")
        for row in selected:
            fwd = object_by_id.get(row["anchor_forward_object_id"])
            rc = object_by_id.get(row["anchor_target_object_id"])
            if fwd is None or rc is None:
                raise ValueError("anchor object absent from object table")

    return object_by_id, relation_by_id, anchors_by_relation, fold_by_cluster


def bounded_credit(residual, target_scale):
    residual = float(residual)
    target_scale = float(target_scale)
    if not np.isfinite(residual) or residual < 0:
        raise ValueError("residual must be finite and non-negative")
    if not np.isfinite(target_scale):
        raise ValueError("target scale must be finite")
    if target_scale <= FORMAL_TARGET_SCALE_FLOOR:
        return None
    return float(target_scale / (target_scale + residual))


def hierarchical_macro(leaves, value_key):
    by_locus = defaultdict(list)
    for row in leaves:
        value = row[value_key]
        if value is not None:
            by_locus[(row["cluster_id"], row["parent_locus_id"])].append(float(value))
    parent_values = {
        (cluster, parent): float(np.mean(values))
        for (cluster, parent), values in by_locus.items()
    }
    by_cluster = defaultdict(list)
    for (cluster, _parent), value in parent_values.items():
        by_cluster[cluster].append(value)
    cluster_values = {cluster: float(np.mean(values)) for cluster, values in by_cluster.items()}
    return (
        float(np.mean(list(cluster_values.values()))) if cluster_values else None,
        parent_values,
        cluster_values,
    )


def bootstrap(cluster_values):
    cluster_ids = sorted(cluster_values)
    if not cluster_ids:
        return None, None
    values = np.asarray([cluster_values[cluster] for cluster in cluster_ids], dtype=np.float64)
    rng = np.random.default_rng(FORMAL_SEED)
    sampled = values[rng.integers(0, len(values), size=(FORMAL_DRAWS, len(values)))].mean(axis=1)
    low, high = np.quantile(sampled, [0.025, 0.975], method="linear")
    return sampled, [float(low), float(high)]


def summarize_leaves(leaves):
    eligible_clusters = {row["cluster_id"] for row in leaves}
    informative = [row for row in leaves if row["F_credit"] is not None]
    informative_clusters = {row["cluster_id"] for row in informative}
    view_fraction = len(informative) / len(leaves) if leaves else 0.0
    cluster_fraction = len(informative_clusters) / len(eligible_clusters) if eligible_clusters else 0.0
    reported = bool(informative_clusters) and view_fraction >= FORMAL_COVERAGE_THRESHOLD and cluster_fraction >= FORMAL_COVERAGE_THRESHOLD
    score, parent_values, cluster_values = hierarchical_macro(informative, "F_credit")
    specificity, _specificity_parents, _specificity_clusters = hierarchical_macro(leaves, "specificity")
    bootstrap_values, interval = bootstrap(cluster_values) if reported else (None, None)
    return {
        "score_status": "REPORTED" if reported else "INSUFFICIENT_DYNAMIC_RANGE",
        "F": score if reported else None,
        "F_95_CI": interval,
        "parent_F": parent_values,
        "cluster_F": cluster_values,
        "F_bootstrap": bootstrap_values,
        "Q_spec": specificity,
        "coverage": {
            "threshold": FORMAL_COVERAGE_THRESHOLD,
            "target_scale_floor": FORMAL_TARGET_SCALE_FLOOR,
            "eligible_views": len(leaves),
            "informative_views": len(informative),
            "view_informative_fraction": view_fraction,
            "eligible_outer_clusters": len(eligible_clusters),
            "informative_outer_clusters": len(informative_clusters),
            "outer_unit_informative_fraction": cluster_fraction,
            "undefined_relation_count": len(leaves) - len(informative),
        },
    }


def publish_result(output, leaves, result, parent_values, cluster_values, bootstrap_values):
    output = Path(output)
    if output.exists():
        raise SystemExit("output path must not already exist")
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        fields = list(leaves[0])
        with (stage / "INDIVIDUAL_SCORES.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(leaves)
        with (stage / "PARENT_F_LEDGER.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["cluster_id", "parent_locus_id", "parent_F"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for (cluster_id, parent_locus_id), value in sorted(parent_values.items()):
                writer.writerow({"cluster_id": cluster_id, "parent_locus_id": parent_locus_id, "parent_F": value})
        with (stage / "DEPENDENCY_CLUSTER_F_LEDGER.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["cluster_id", "cluster_F"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for cluster_id, value in sorted(cluster_values.items()):
                writer.writerow({"cluster_id": cluster_id, "cluster_F": value})
        if bootstrap_values is not None:
            np.save(stage / "F_BOOTSTRAP.npy", np.asarray(bootstrap_values, dtype=np.float64), allow_pickle=False)
        (stage / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        names = ["INDIVIDUAL_SCORES.tsv", "PARENT_F_LEDGER.tsv", "DEPENDENCY_CLUSTER_F_LEDGER.tsv", "FORMAL_RESULT.json"]
        if bootstrap_values is not None:
            names.append("F_BOOTSTRAP.npy")
        manifest = {name: {"sha256": sha256_file(stage / name), "bytes": (stage / name).stat().st_size} for name in names}
        (stage / "OUTPUT_MANIFEST.json").write_text(json.dumps({"status": "PASS", "files": manifest}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        names.append("OUTPUT_MANIFEST.json")
        (stage / "SHA256SUMS").write_text("".join(f"{sha256_file(stage / name)}  {name}\n" for name in names), encoding="utf-8")
        (stage / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")
        os.replace(stage, output)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
