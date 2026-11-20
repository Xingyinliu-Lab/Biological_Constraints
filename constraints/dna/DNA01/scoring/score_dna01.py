#!/usr/bin/env python3
"""DNA01 v0.5.4 public F-only full-chain shared-affine scorer."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from dna01_core import (
    FORMAL_ALPHA,
    FORMAL_DRAWS,
    FORMAL_SEED,
    bounded_credit,
    load_vectors,
    publish_result,
    read_tsv,
    sha256_file,
    summarize_leaves,
    validate_refset,
)


DNA01_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFSET = DNA01_ROOT / "refset"


def hierarchical_weights(train):
    by_cluster = defaultdict(list)
    for row in train:
        by_cluster[row["source_dependency_cluster_id"]].append(row)
    weights = {}
    for cluster_rows in by_cluster.values():
        by_parent = defaultdict(list)
        for row in cluster_rows:
            by_parent[row["parent_locus_id"]].append(row)
        for parent_rows in by_parent.values():
            weight = 1.0 / len(by_cluster) / len(by_parent) / len(parent_rows)
            for row in parent_rows:
                weights[row["relation_id"]] = weight
    return weights


def fit_affine(train, embeddings):
    weights = hierarchical_weights(train)
    x = np.stack([embeddings[row["forward_object_id"]] for row in train])
    y = np.stack([embeddings[row["rc_object_id"]] for row in train])
    w = np.asarray([weights[row["relation_id"]] for row in train], dtype=np.float64)
    mean_x = np.sum(x * w[:, None], axis=0)
    mean_y = np.sum(y * w[:, None], axis=0)
    centered_x = x - mean_x
    centered_y = y - mean_y
    dimension = x.shape[1]
    scale = float(np.sum(w * np.sum(centered_x * centered_x, axis=1)) / dimension)
    ridge = FORMAL_ALPHA * max(scale, np.finfo(np.float64).eps)
    gram = (centered_x.T * w) @ centered_x + ridge * np.eye(dimension)
    rhs = (centered_x.T * w) @ centered_y
    matrix = np.linalg.solve(gram, rhs)
    offset = mean_y - mean_x @ matrix
    if not np.isfinite(matrix).all() or not np.isfinite(offset).all():
        raise ValueError("affine fit produced non-finite parameters")
    return matrix, offset, ridge


def main():
    parser = argparse.ArgumentParser(
        description="Score a model representation against DNA01.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--relations", type=Path, default=DEFAULT_REFSET / "DNA01_RELATION_TABLE.tsv")
    parser.add_argument("--anchors", type=Path, default=DEFAULT_REFSET / "DNA01_ANCHOR_TABLE.tsv")
    parser.add_argument("--objects", type=Path, default=DEFAULT_REFSET / "DNA01_OBJECT_TABLE.tsv")
    parser.add_argument("--fold-assignment", type=Path, default=DEFAULT_REFSET / "DNA01_FOLD_ASSIGNMENT_TABLE.tsv")
    parser.add_argument("--dependency-clusters", type=Path, default=DEFAULT_REFSET / "DNA01_DEPENDENCY_CLUSTER_TABLE.tsv")
    parser.add_argument("--repeat-control-panel", type=Path, default=DEFAULT_REFSET / "DNA01_REPEAT_CONTROL_PANEL.tsv")
    parser.add_argument("--embeddings", required=True, type=Path)
    parser.add_argument("--repeat-embeddings", required=True, type=Path)
    parser.add_argument("--model-key", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    relations = read_tsv(args.relations)
    anchors = read_tsv(args.anchors)
    objects = read_tsv(args.objects)
    folds = read_tsv(args.fold_assignment)
    dependency_rows = read_tsv(args.dependency_clusters)
    repeat_rows = read_tsv(args.repeat_control_panel)
    object_by_id, _relation_by_id, anchors_by_relation, _fold_by_cluster = validate_refset(
        objects, relations, anchors, folds, dependency_rows, repeat_rows
    )

    primary_raw, primary_dtype = load_vectors(args.embeddings)
    repeat_raw, repeat_dtype = load_vectors(args.repeat_embeddings)
    formal_ids = set(object_by_id)
    if set(primary_raw) != formal_ids or set(repeat_raw) != formal_ids:
        raise ValueError("primary and repeat embeddings must exactly cover the formal object universe")
    if next(iter(primary_raw.values())).shape != next(iter(repeat_raw.values())).shape:
        raise ValueError("primary and repeat representation dimensions differ")
    if primary_dtype != repeat_dtype:
        raise ValueError("primary and repeat representation dtypes differ")
    mismatched = [
        oid for oid in sorted(formal_ids)
        if primary_raw[oid].tobytes(order="C") != repeat_raw[oid].tobytes(order="C")
    ]
    if mismatched:
        raise ValueError(f"repeat forward is not bitwise deterministic for {len(mismatched)} objects")
    embeddings = {oid: np.asarray(vector, dtype=np.float64) for oid, vector in primary_raw.items()}

    folds_present = sorted({int(row["outer_fold_id"]) for row in relations})
    if folds_present != [0, 1, 2, 3, 4]:
        raise ValueError("expected five outer folds")
    operators = {}
    ridge_by_fold = {}
    for fold in folds_present:
        training = [row for row in relations if int(row["outer_fold_id"]) != fold]
        operators[fold] = fit_affine(training, embeddings)
        ridge_by_fold[str(fold)] = operators[fold][2]

    leaves = []
    for relation in relations:
        fold = int(relation["outer_fold_id"])
        matrix, offset, _ridge = operators[fold]
        forward = embeddings[relation["forward_object_id"]]
        target = embeddings[relation["rc_object_id"]]
        predicted = forward @ matrix + offset
        selected = anchors_by_relation[relation["relation_id"]]
        wrong_distances = [float(np.linalg.norm(predicted - embeddings[row["anchor_target_object_id"]])) for row in selected]
        target_distances = [float(np.linalg.norm(target - embeddings[row["anchor_target_object_id"]])) for row in selected]
        target_scale = float(np.median(target_distances))
        residual = float(np.linalg.norm(predicted - target))
        credit = bounded_credit(residual, target_scale)
        specificity = float(np.mean([1.0 if residual < value else 0.5 if residual == value else 0.0 for value in wrong_distances]))
        cycle_residual = float(np.linalg.norm((predicted @ matrix + offset) - forward))
        cycle_scale = float(np.median([np.linalg.norm(forward - embeddings[row["anchor_forward_object_id"]]) for row in selected]))
        normalized_cycle = cycle_residual / cycle_scale if cycle_scale > 0 else None
        leaf = {
            "relation_id": relation["relation_id"],
            "forward_object_id": relation["forward_object_id"],
            "rc_object_id": relation["rc_object_id"],
            "outer_fold_id": fold,
            "cluster_id": relation["source_dependency_cluster_id"],
            "parent_locus_id": relation["parent_locus_id"],
            "view_id": relation["view_id"],
            "evaluation_status": "SCORED" if credit is not None else "UNSCORABLE",
            "unscorable_reason": "" if credit is not None else "ZERO_TARGET_SCALE",
            "residual": residual,
            "target_scale": target_scale,
            "F_credit": credit,
            "specificity": specificity,
            "cycle_residual": cycle_residual,
            "cycle_scale": cycle_scale,
            "normalized_cycle": normalized_cycle,
        }
        for rank, distance in enumerate(wrong_distances, start=1):
            leaf[f"anchor_target_distance_{rank}"] = distance
        for rank, distance in enumerate(target_distances, start=1):
            leaf[f"target_to_anchor_distance_{rank}"] = distance
        leaves.append(leaf)

    summary = summarize_leaves(leaves)
    cycles = [row["normalized_cycle"] for row in leaves if row["normalized_cycle"] is not None]
    result = {
        "status": "PASS",
        "score_status": summary["score_status"],
        "constraint_id": "C-DNA-01",
        "model_key": args.model_key,
        "result_kind": "MODEL_REPRESENTATION",
        "operator_branch": "EMPIRICAL_SHARED_AFFINE",
        "score_transform_version": "DNA01_v0.5.4_F_ONLY_FULL_CHAIN",
        "score_transform": "F_credit = a_i/(a_i+e_i)",
        "F": summary["F"],
        "F_95_CI": summary["F_95_CI"],
        "Q_spec": summary["Q_spec"],
        "representation_status": "EMPIRICAL_AFFINE_CHARACTERIZATION",
        "equivariance_validity": "UNDETERMINED",
        "numeric_method": {
            "alpha": FORMAL_ALPHA,
            "alpha_scaling": "alpha * weighted_centered_input_energy_per_dimension",
            "outer_folds": 5,
            "fit_population": "four outer folds; held-out fold excluded",
            "bootstrap": {"draws": FORMAL_DRAWS, "seed": FORMAL_SEED, "unit": "dependency_cluster", "operator_refit": False},
            "execution_mode": "FORMAL",
        },
        "coverage": summary["coverage"],
        "counts": {"objects": len(objects), "relations": len(relations), "anchors": len(anchors), "clusters": len(folds)},
        "repeat_forward": {"bitwise_identical": True, "nonidentical_object_count": 0},
        "cycle": {
            "median_normalized_residual": float(np.median(cycles)) if cycles else None,
            "maximum_normalized_residual": float(np.max(cycles)) if cycles else None,
        },
        "ridge_by_outer_fold": ridge_by_fold,
        "input_sha256": {
            "relations": sha256_file(args.relations),
            "anchors": sha256_file(args.anchors),
            "objects": sha256_file(args.objects),
            "fold_assignment": sha256_file(args.fold_assignment),
            "dependency_clusters": sha256_file(args.dependency_clusters),
            "repeat_control_panel": sha256_file(args.repeat_control_panel),
            "embeddings": sha256_file(args.embeddings),
            "repeat_embeddings": sha256_file(args.repeat_embeddings),
        },
    }
    publish_result(args.output, leaves, result, summary["parent_F"], summary["cluster_F"], summary["F_bootstrap"])
    print(json.dumps({"model": args.model_key, "F": result["F"], "F_95_CI": result["F_95_CI"], "Q_spec": result["Q_spec"]}, sort_keys=True))


if __name__ == "__main__":
    main()
