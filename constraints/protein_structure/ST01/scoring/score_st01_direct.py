"""Fail-closed direct intrinsic-ranking scorer for C-ST-01."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_relations(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"contrast_id", "positive_edge_id", "anchor_object_id", "true_target_object_id", "negative_object_id", "source_dependency_cluster_id"}
    if not rows or set(rows[0]) < required:
        raise ValueError("relations table lacks required fields")
    if len({row["contrast_id"] for row in rows}) != len(rows):
        raise ValueError("duplicate contrast_id")
    edge_anchors: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        edge_anchors[row["positive_edge_id"]].add(row["anchor_object_id"])
    if any(len(anchors) != 2 for anchors in edge_anchors.values()) or any(sum(row["positive_edge_id"] == edge_id for row in rows) != 2 for edge_id in edge_anchors):
        raise ValueError("strict edge bundle is incomplete")
    return rows


def aggregate(rows: list[dict[str, str]], pair_score: dict[tuple[str, str], float], tie_tolerance: float, seed: int, draws: int) -> tuple[dict, dict[str, list[dict[str, object]]], np.ndarray]:
    endpoint_credit: dict[tuple[str, str], list[float]] = defaultdict(list)
    endpoint_cluster: dict[tuple[str, str], str] = {}
    contrast_rows: list[dict[str, object]] = []
    for row in rows:
        anchor, true_target, negative = row["anchor_object_id"], row["true_target_object_id"], row["negative_object_id"]
        try:
            true_score = pair_score[(anchor, true_target)]
            negative_score = pair_score[(anchor, negative)]
            delta = true_score - negative_score
        except KeyError as error:
            raise ValueError("missing required pair score") from error
        if not np.isfinite(delta):
            raise ValueError("nonfinite formal pair score")
        f_credit = 1.0 if delta > tie_tolerance else -1.0 if delta < -tie_tolerance else 0.0
        contrast_rows.append({
            "contrast_id": row["contrast_id"],
            "positive_edge_id": row["positive_edge_id"],
            "anchor_object_id": anchor,
            "true_target_object_id": true_target,
            "negative_object_id": negative,
            "source_dependency_cluster_id": row["source_dependency_cluster_id"],
            "true_pair_score": true_score,
            "negative_pair_score": negative_score,
            "score_delta": delta,
            "F_credit": f_credit,
        })
        endpoint = (row["positive_edge_id"], anchor)
        endpoint_credit[endpoint].append(f_credit)
        prior = endpoint_cluster.setdefault(endpoint, row["source_dependency_cluster_id"])
        if prior != row["source_dependency_cluster_id"]:
            raise ValueError("endpoint maps to multiple dependency clusters")

    endpoint_rows: list[dict[str, object]] = []
    anchor_credit: dict[tuple[str, str], list[float]] = defaultdict(list)
    for endpoint, credits in endpoint_credit.items():
        value = float(np.mean(credits))
        endpoint_rows.append({"positive_edge_id": endpoint[0], "anchor_object_id": endpoint[1], "source_dependency_cluster_id": endpoint_cluster[endpoint], "endpoint_F": value, "contrast_count": len(credits)})
        anchor_credit[(endpoint[1], endpoint_cluster[endpoint])].append(value)
    anchor_rows: list[dict[str, object]] = []
    cluster_credit: dict[str, list[float]] = defaultdict(list)
    for (anchor, cluster), credits in anchor_credit.items():
        value = float(np.mean(credits))
        anchor_rows.append({"anchor_object_id": anchor, "source_dependency_cluster_id": cluster, "anchor_F": value, "endpoint_count": len(credits)})
        cluster_credit[cluster].append(value)
    cluster_rows = [{"source_dependency_cluster_id": cluster, "cluster_F": float(np.mean(cluster_credit[cluster])), "anchor_count": len(cluster_credit[cluster])} for cluster in sorted(cluster_credit)]
    cluster_values = np.array([np.mean(cluster_credit[cluster]) for cluster in sorted(cluster_credit)], dtype=float)
    if not np.all(np.isfinite(cluster_values)) or len(cluster_values) == 0:
        raise ValueError("invalid cluster aggregation")
    f_value = float(np.mean(cluster_values))
    rng = np.random.Generator(np.random.PCG64(seed))
    bootstrap = cluster_values[rng.integers(0, len(cluster_values), size=(draws, len(cluster_values)))].mean(axis=1)
    ci = np.quantile(bootstrap, [0.025, 0.975], method="linear")
    result = {
        "formal_interface": "F_ONLY",
        "score_transform_version": "ST01_v0.1.2_F_ONLY",
        "F": f_value,
        "F_95_CI": [float(ci[0]), float(ci[1])],
        "counts": {"directed_contrasts": len(rows), "endpoint_challenges": len(endpoint_credit), "anchors": len(anchor_credit), "source_dependency_clusters": len(cluster_values)},
    }
    ledgers = {"contrast": contrast_rows, "endpoint": endpoint_rows, "anchor": anchor_rows, "cluster": cluster_rows}
    return result, ledgers, bootstrap


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty ledger: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def embedding_pair_scores(path: Path, rows: list[dict[str, str]]) -> dict[tuple[str, str], float]:
    data = np.load(path, allow_pickle=False)
    if set(data.files) != {"object_ids", "embeddings"}:
        raise ValueError("embedding NPZ must contain exactly object_ids and embeddings")
    object_ids = [str(x) for x in data["object_ids"]]
    embeddings = np.asarray(data["embeddings"], dtype=float)
    if embeddings.ndim != 2 or embeddings.shape[1] <= 0 or embeddings.shape[0] != len(object_ids) or len(set(object_ids)) != len(object_ids):
        raise ValueError("invalid embedding matrix or IDs")
    norm = np.linalg.norm(embeddings, axis=1)
    if not np.all(np.isfinite(embeddings)) or np.any(norm == 0):
        raise ValueError("nonfinite or zero-norm embedding")
    normalized = embeddings / norm[:, None]
    index = {object_id: i for i, object_id in enumerate(object_ids)}
    needed = {row[key] for row in rows for key in ("anchor_object_id", "true_target_object_id", "negative_object_id")}
    if not needed <= set(index):
        raise ValueError("embedding artifact misses formal object")
    return {(anchor, candidate): float(np.dot(normalized[index[anchor]], normalized[index[candidate]])) for anchor in needed for candidate in needed if anchor != candidate}


def ledger_pair_scores(path: Path, rows: list[dict[str, str]]) -> dict[tuple[str, str], float]:
    with path.open(encoding="utf-8", newline="") as handle:
        ledger = list(csv.DictReader(handle, delimiter="\t"))
    scores = {}
    for row in ledger:
        key = (row["anchor_object_id"], row["candidate_object_id"])
        if key in scores:
            raise ValueError("duplicate comparator pair score")
        score = float(row["score"])
        if not np.isfinite(score):
            raise ValueError("nonfinite comparator pair score")
        scores[key] = score
    return scores


def main() -> None:
    constraint_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", type=Path, default=constraint_root / "refset" / "relations.tsv")
    parser.add_argument("--numeric-method", type=Path, default=constraint_root / "scoring" / "ST01_NUMERIC_METHOD_v0.1.1_FROZEN_20261010.json")
    parser.add_argument("--binding", type=Path, default=constraint_root / "scoring" / "ST01_NUMERIC_INPUT_BINDING.json")
    parser.add_argument("--embeddings", type=Path)
    parser.add_argument("--pair-scores", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for the complete public output package")
    args = parser.parse_args()
    if bool(args.embeddings) == bool(args.pair_scores):
        raise ValueError("provide exactly one of --embeddings or --pair-scores")
    method, binding = json.loads(args.numeric_method.read_text(encoding="utf-8")), json.loads(args.binding.read_text(encoding="utf-8"))
    if method["primary_evaluation"]["train_test_split"] != "none" or method["primary_evaluation"]["trained_parameters"] is not False:
        raise ValueError("Numeric Method is incompatible with direct ST01 scoring")
    if binding["numeric_method_sha256"] != sha256_file(args.numeric_method):
        raise ValueError("numeric method hash does not match input binding")
    if binding["relation_table_sha256"] != sha256_file(args.relations):
        raise ValueError("relation table hash does not match input binding")
    rows = load_relations(args.relations)
    pair_scores = embedding_pair_scores(args.embeddings, rows) if args.embeddings else ledger_pair_scores(args.pair_scores, rows)
    result, ledgers, bootstrap = aggregate(rows, pair_scores, float(method["tie_tolerance"]), int(method["bootstrap"]["seed"]), int(method["bootstrap"]["draws"]))
    forward_input = args.embeddings if args.embeddings else args.pair_scores
    result.update({
        "numeric_method_sha256": sha256_file(args.numeric_method),
        "input_binding_sha256": sha256_file(args.binding),
        "relation_table_sha256": sha256_file(args.relations),
        "scorer_id": "ST01_DIRECT_INTRINSIC_RANKING_SCORER_v0.1.1",
        "forward_input": {
            "mode": "embeddings_npz" if args.embeddings else "pair_scores_tsv",
            "path": str(forward_input),
            "sha256": sha256_file(forward_input),
        },
    })
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(output_dir / "INDIVIDUAL_SCORES.tsv", ledgers["contrast"])
    write_tsv(output_dir / "ENDPOINT_SCORES.tsv", ledgers["endpoint"])
    write_tsv(output_dir / "ANCHOR_SCORES.tsv", ledgers["anchor"])
    write_tsv(output_dir / "CLUSTER_SCORES.tsv", ledgers["cluster"])
    np.save(output_dir / "F_BOOTSTRAP.npy", bootstrap, allow_pickle=False)
    (output_dir / "SCORE_SUMMARY.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {}
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "OUTPUT_MANIFEST.json":
            manifest[path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (output_dir / "OUTPUT_MANIFEST.json").write_text(json.dumps({"formal_interface": "F_ONLY", "score_transform_version": "ST01_v0.1.2_F_ONLY", "files": manifest}, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
