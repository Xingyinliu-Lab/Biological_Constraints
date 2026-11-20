from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from rna04_public_core import aggregate, site_scores, validate_refset

BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 20261002
C_CANDIDATES = (1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_vectors(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {"sequence_sha256", "vectors"}:
            raise ValueError("vector NPZ must contain exactly sequence_sha256 and vectors")
        raw_ids = data["sequence_sha256"]
        raw_values = data["vectors"]
        if raw_ids.ndim != 1 or raw_ids.dtype.kind not in "US":
            raise ValueError("sequence_sha256 must be a one-dimensional string array")
        if raw_values.dtype.kind not in "iuf" or np.issubdtype(raw_values.dtype, np.bool_):
            raise ValueError("vectors must contain real numeric values; complex, boolean, string and object arrays are invalid")
        ids = [value.decode() if isinstance(value, bytes) else str(value) for value in raw_ids]
        values = np.asarray(raw_values, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] <= 0 or values.shape[0] != len(ids) or len(set(ids)) != len(ids):
        raise ValueError("invalid vector matrix or duplicate sequence identity")
    if not np.isfinite(values).all():
        raise ValueError("non-finite vector")
    return dict(zip(ids, values))


def training_weights(rows: list[dict]) -> np.ndarray:
    hierarchy = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    for index, row in enumerate(rows):
        hierarchy[row["outer_cluster_id"]][row["parent_object_id"]][row["site_id"]][row["mutation_class"]].append(index)
    result = np.zeros(len(rows))
    cluster_count = len(hierarchy)
    for parents in hierarchy.values():
        for sites in parents.values():
            for classes in sites.values():
                for indices in classes.values():
                    result[indices] = 1 / cluster_count / len(parents) / len(sites) / len(classes) / len(indices)
    return result


def group_folds(groups: np.ndarray, count: int) -> list[tuple[np.ndarray, np.ndarray]]:
    unique = sorted(set(groups))
    order = sorted(unique, key=lambda value: hashlib.sha256(value.encode()).hexdigest())
    folds: list[list[str]] = [[] for _ in range(count)]
    for index, group in enumerate(order):
        folds[index % count].append(group)
    result = []
    for values in folds:
        selected = set(values)
        test = np.asarray([index for index, group in enumerate(groups) if group in selected], dtype=int)
        train = np.asarray([index for index, group in enumerate(groups) if group not in selected], dtype=int)
        result.append((train, test))
    return result


def score_split(train: np.ndarray, test: np.ndarray, features: np.ndarray, rows: list[dict], regularization: float) -> np.ndarray:
    mean = features[train].mean(axis=0)
    scale = features[train].std(axis=0)
    scale[scale == 0] = 1.0
    train_x = (features[train] - mean) / scale
    test_x = (features[test] - mean) / scale
    labels = np.asarray([1.0 if rows[index]["mutation_class"] == "missense" else 0.0 for index in train])
    weights = training_weights([rows[index] for index in train])
    weights = weights / weights.sum()
    alpha = 1.0 / regularization

    def objective(theta: np.ndarray):
        coefficients, intercept = theta[:-1], theta[-1]
        linear = train_x @ coefficients + intercept
        loss = np.sum(weights * (np.logaddexp(0.0, linear) - labels * linear)) + 0.5 * alpha * np.dot(coefficients, coefficients)
        probability = 1.0 / (1.0 + np.exp(-np.clip(linear, -40, 40)))
        derivative = (probability - labels) * weights
        return loss, np.r_[train_x.T @ derivative + alpha * coefficients, derivative.sum()]

    fitted = minimize(lambda theta: objective(theta), np.zeros(train_x.shape[1] + 1), jac=True, method="L-BFGS-B", options={"maxiter": 300, "ftol": 1e-9})
    if not fitted.success:
        raise RuntimeError(f"ridge optimizer failed: {fitted.message}")
    linear = test_x @ fitted.x[:-1] + fitted.x[-1]
    return 1.0 / (1.0 + np.exp(-np.clip(linear, -40, 40)))


def main() -> None:
    parser = argparse.ArgumentParser(description="RNA04 public nested-cross-fit scorer")
    refset = Path(__file__).resolve().parents[1] / "refset"
    default_ledger = refset / "COMPLETE_SITE_ALT_LEDGER.tsv"
    parser.add_argument("--ledger", type=Path, default=default_ledger)
    parser.add_argument("--fasta", type=Path, default=refset / "RNA04_unique_endpoints.fasta")
    parser.add_argument("--vectors", type=Path, required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = read_tsv(args.ledger)
    refset_identity = validate_refset(rows, args.fasta)
    registered_identity = validate_refset(read_tsv(default_ledger), refset / "RNA04_unique_endpoints.fasta")["population_identity_sha256"]
    vectors = load_vectors(args.vectors)
    required = {row[key] for row in rows for key in ("wt_sequence_sha256", "alt_sequence_sha256")}
    if set(vectors) != required:
        raise ValueError(f"vector coverage mismatch; missing={len(required-set(vectors))}, extra={len(set(vectors)-required)}")
    features = np.stack([vectors[row["alt_sequence_sha256"]] - vectors[row["wt_sequence_sha256"]] for row in rows])
    if not np.isfinite(features).all():
        raise ValueError("non-finite feature matrix")
    groups = np.asarray([row["outer_cluster_id"] for row in rows])
    predictions = np.full(len(rows), np.nan)
    selected_parameters = []
    for fold, (train, test) in enumerate(group_folds(groups, 5), 1):
        training_groups = groups[train]
        best = None
        for regularization in C_CANDIDATES:
            validation_scores = []
            for inner_train, inner_validation in group_folds(training_groups, 4):
                inner_train_indices = train[inner_train]
                validation_indices = train[inner_validation]
                scores = score_split(inner_train_indices, validation_indices, features, rows, regularization)
                sites, _, _ = site_scores([rows[index] for index in validation_indices], scores)
                validation_scores.append(aggregate([rows[index] for index in validation_indices], sites)[0])
            key = (float(np.mean(validation_scores)), -regularization)
            if best is None or key > best[0]:
                best = (key, regularization)
        regularization = best[1]
        predictions[test] = score_split(train, test, features, rows, regularization)
        selected_parameters.append({
            "outer_fold": fold, "C": regularization, "inner_hierarchical_F": best[0][0],
            "test_rows": len(test), "test_clusters": len(set(groups[test])),
        })
    if not np.isfinite(predictions).all():
        raise ValueError("missing held-out scores")

    site_values, site_records, pair_records = site_scores(rows, predictions)
    formal_F, cluster_values, parent_values = aggregate(rows, site_values)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    cluster_ids = np.asarray(sorted(cluster_values))
    cluster_array = np.asarray([cluster_values[key] for key in cluster_ids])
    bootstrap_means = np.mean(rng.choice(cluster_array, size=(BOOTSTRAP_DRAWS, len(cluster_array)), replace=True), axis=1)
    low, high = np.quantile(bootstrap_means, [0.025, 0.975])

    args.output.mkdir(parents=True, exist_ok=True)
    prediction_rows = [{
        "row_id": index, "site_id": row["site_id"], "parent_object_id": row["parent_object_id"],
        "outer_cluster_id": row["outer_cluster_id"], "mutation_class": row["mutation_class"],
        "wt_sequence_sha256": row["wt_sequence_sha256"], "alt_sequence_sha256": row["alt_sequence_sha256"],
        "heldout_probability_missense": format(float(predictions[index]), ".17g"),
    } for index, row in enumerate(rows)]
    parent_rows = [{"parent_object_id": key, "parent_F": format(parent_values[key], ".17g")} for key in sorted(parent_values)]
    cluster_rows = [{"outer_cluster_id": key, "cluster_F": format(cluster_values[key], ".17g")} for key in sorted(cluster_values)]
    paths = {
        "HELDOUT_PREDICTIONS.tsv": prediction_rows,
        "PAIR_CREDIT_LEDGER.tsv": pair_records,
        "INDIVIDUAL_SCORES.tsv": site_records,
        "PARENT_SCORE_LEDGER.tsv": parent_rows,
        "CLUSTER_SCORE_LEDGER.tsv": cluster_rows,
    }
    for name, records in paths.items():
        write_tsv(args.output / name, records)
    result_path = args.output / "FORMAL_RESULT.json"
    result_path.write_text(json.dumps({
        "status": "PASS", "formal_result": True, "assay_id": "C-RNA-04", "model_id": args.model_id,
        "feature": "whole-CDS vector delta (ALT minus WT)",
        "formal_interface": "F_ONLY", "score_transform_version": "RNA04_v0.2.2_F_ONLY",
        "readout": "outer-CV L2 logistic ridge; inner hierarchical site-parent-cluster F selection",
        "training_weights": "equal outer cluster -> parent -> site -> class -> variant", "class_weight": "none",
        "rows": len(rows), "sites": len(site_values), "parents": len(parent_values), "outer_clusters": len(cluster_values),
        "F": formal_F, "F_95_CI": [float(low), float(high)],
        "bootstrap_replicates": BOOTSTRAP_DRAWS, "bootstrap_seed": BOOTSTRAP_SEED, "bootstrap_unit": "outer_cluster",
        "outer_fold_C": selected_parameters, "ledger_sha256": sha256(args.ledger), "fasta_sha256": sha256(args.fasta),
        "population_identity_sha256": refset_identity["population_identity_sha256"],
        "registered_reference_population": refset_identity["population_identity_sha256"] == registered_identity,
        "vectors_sha256": sha256(args.vectors),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    np.save(args.output / "BOOTSTRAP_MEANS.npy", bootstrap_means)
    output_files = [args.output / name for name in paths] + [result_path, args.output / "BOOTSTRAP_MEANS.npy"]
    (args.output / "OUTPUT_MANIFEST.json").write_text(json.dumps({
        "files": [{"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in output_files]
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
