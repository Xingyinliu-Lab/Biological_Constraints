"""Run the public frozen ST02 scoring protocol on structure representations."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from st02_numeric_core import aggregate_hierarchy, bayesian_cluster_bootstrap

REQUIRED_RELATION_FIELDS = {
    "contrast_id", "positive_edge_id", "direction_id", "anchor_object_id",
    "true_target_object_id", "negative_object_id", "topology_id",
    "source_homologous_superfamily_id", "target_homologous_superfamily_id",
    "source_dependency_cluster_id",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_relations(path: Path, expected_k: int) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not REQUIRED_RELATION_FIELDS.issubset(rows[0]):
        raise ValueError("relation schema failure")
    contrast_ids = [row["contrast_id"] for row in rows]
    if len(contrast_ids) != len(set(contrast_ids)):
        raise ValueError("duplicate contrast_id")
    directions: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        directions.setdefault(row["direction_id"], []).append(row)
    if any(len(group) != expected_k for group in directions.values()):
        raise ValueError("directed endpoint K closure failure")
    if any(len({row["negative_object_id"] for row in group}) != expected_k for group in directions.values()):
        raise ValueError("duplicate negative within endpoint")
    edge_directions: dict[str, set[str]] = {}
    for row in rows:
        edge_directions.setdefault(row["positive_edge_id"], set()).add(row["direction_id"])
    if any(len(direction_ids) != 2 for direction_ids in edge_directions.values()):
        raise ValueError("bidirectional edge closure failure")
    return rows


def stringify(mapping: dict) -> dict:
    return {"|".join(key) if isinstance(key, tuple) else key: value for key, value in mapping.items()}


def mean(values: list[float]) -> float:
    return float(np.mean(np.asarray(values, dtype=np.float64)))


def hierarchy_ledgers(rows: list[dict[str, str]], credits: list[float]) -> dict[str, list[dict[str, object]]]:
    endpoint_values: dict[str, list[float]] = defaultdict(list)
    endpoint_meta: dict[str, tuple[str, str, str, str, str]] = {}
    for row, credit in zip(rows, credits, strict=True):
        direction = row["direction_id"]
        endpoint_values[direction].append(float(credit))
        endpoint_meta[direction] = (
            row["anchor_object_id"], row["target_homologous_superfamily_id"],
            row["source_dependency_cluster_id"], row["source_homologous_superfamily_id"], row["topology_id"],
        )
    endpoint_rows = []
    target_values: dict[tuple[str, str, str, str, str], list[float]] = defaultdict(list)
    for direction in sorted(endpoint_values):
        value = mean(endpoint_values[direction])
        anchor, target_sf, cluster, source_sf, topology = endpoint_meta[direction]
        endpoint_rows.append({"direction_id": direction, "anchor_object_id": anchor, "target_homologous_superfamily_id": target_sf, "source_dependency_cluster_id": cluster, "source_homologous_superfamily_id": source_sf, "topology_id": topology, "endpoint_F": value, "contrast_count": len(endpoint_values[direction])})
        target_values[(anchor, target_sf, cluster, source_sf, topology)].append(value)
    target_rows = []
    anchor_values: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    for key in sorted(target_values):
        anchor, target_sf, cluster, source_sf, topology = key
        value = mean(target_values[key])
        target_rows.append({"anchor_object_id": anchor, "target_homologous_superfamily_id": target_sf, "source_dependency_cluster_id": cluster, "source_homologous_superfamily_id": source_sf, "topology_id": topology, "target_sf_F": value, "endpoint_count": len(target_values[key])})
        anchor_values[(anchor, cluster, source_sf, topology)].append(value)
    anchor_rows = []
    cluster_values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for key in sorted(anchor_values):
        anchor, cluster, source_sf, topology = key
        value = mean(anchor_values[key])
        anchor_rows.append({"anchor_object_id": anchor, "source_dependency_cluster_id": cluster, "source_homologous_superfamily_id": source_sf, "topology_id": topology, "anchor_F": value, "target_sf_count": len(anchor_values[key])})
        cluster_values[(cluster, source_sf, topology)].append(value)
    cluster_rows = []
    for key in sorted(cluster_values):
        cluster, source_sf, topology = key
        cluster_rows.append({"source_dependency_cluster_id": cluster, "source_homologous_superfamily_id": source_sf, "topology_id": topology, "cluster_F": mean(cluster_values[key]), "anchor_count": len(cluster_values[key])})
    return {"endpoint": endpoint_rows, "target_sf": target_rows, "anchor": anchor_rows, "cluster": cluster_rows}


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"empty public ledger: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    constraint_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--refset-dir", type=Path, default=constraint_root / "refset")
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="Complete formal F-only output directory")
    args = parser.parse_args()
    output_path = args.output_dir / "SCORE_SUMMARY.json"
    if output_path.exists():
        raise FileExistsError(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    method_path = Path(__file__).with_name("ST02_NUMERIC_METHOD_v0.1.3_F_ONLY.json")
    method = json.loads(method_path.read_text(encoding="utf-8"))
    if method.get("status") != "FROZEN_NUMERIC_METHOD":
        raise ValueError("Numeric Method is not frozen")
    uncertainty = method["uncertainty"]
    if uncertainty["scheme"] != "GLOBAL_DEPENDENCY_CLUSTER_BAYESIAN_BOOTSTRAP":
        raise ValueError("unsupported bootstrap scheme")

    relations_path = args.refset_dir / "relations.tsv"
    rows = read_relations(relations_path, expected_k=4)
    needed = {
        row[field]
        for row in rows
        for field in ("anchor_object_id", "true_target_object_id", "negative_object_id")
    }
    with np.load(args.embeddings, allow_pickle=False) as archive:
        if set(archive.files) != {"object_ids", "embeddings"}:
            raise ValueError("embedding NPZ schema failure")
        object_ids = [str(value) for value in archive["object_ids"].tolist()]
        matrix = np.asarray(archive["embeddings"], dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[1] <= 0 or matrix.shape[0] != len(object_ids) or len(object_ids) != len(set(object_ids)):
        raise ValueError("embedding shape/identity failure")
    if set(object_ids) != needed:
        raise ValueError("embedding population mismatch")
    if not np.isfinite(matrix).all():
        raise ValueError("nonfinite embedding")
    norms = np.linalg.norm(matrix, axis=1)
    if np.any(norms == 0):
        raise ValueError("zero-norm embedding")
    matrix /= norms[:, None]
    index = {object_id: row_index for row_index, object_id in enumerate(object_ids)}

    credits = []
    individual_rows = []
    for row in rows:
        anchor = matrix[index[row["anchor_object_id"]]]
        positive = float(anchor @ matrix[index[row["true_target_object_id"]]])
        negative = float(anchor @ matrix[index[row["negative_object_id"]]])
        delta = positive - negative
        credit = 1.0 if delta > 0 else -1.0 if delta < 0 else 0.0
        credits.append(credit)
        individual_rows.append({
            "contrast_id": row["contrast_id"], "positive_edge_id": row["positive_edge_id"], "direction_id": row["direction_id"],
            "anchor_object_id": row["anchor_object_id"], "true_target_object_id": row["true_target_object_id"], "negative_object_id": row["negative_object_id"],
            "topology_id": row["topology_id"], "source_homologous_superfamily_id": row["source_homologous_superfamily_id"],
            "target_homologous_superfamily_id": row["target_homologous_superfamily_id"], "source_dependency_cluster_id": row["source_dependency_cluster_id"],
            "positive_cosine": positive, "negative_cosine": negative, "score_delta": delta, "F_credit": credit,
        })

    point = aggregate_hierarchy(rows, credits)
    bootstrap = bayesian_cluster_bootstrap(
        rows,
        credits,
        draws=int(uncertainty["draws"]),
        seed=int(uncertainty["seed"]),
        return_overall_draws=True,
    )
    overall_draws = bootstrap.pop("overall_draws", None)
    result = {
        "artifact_id": "ST02_F_ONLY_FORMAL_EVALUATION",
        "status": "PASS",
        "evaluation_class": "FORMAL_PROTOCOL_EVALUATION",
        "formal_scoring_protocol": True,
        "registered_paper_result": False,
        "model_id": args.model_id,
        "formal_interface": "F_ONLY",
        "score_transform_version": "ST02_v0.1.3_F_ONLY",
        "F": point["F"],
        "F_95_CI": bootstrap["overall_ci95"],
        "topology_F": stringify(point["topology"]),
        "topology_F_95_CI": stringify(bootstrap["topology_ci95"]),
        "source_sf_topology_F": stringify(point["source_sf_topology"]),
        "source_sf_topology_F_95_CI": stringify(bootstrap["source_sf_topology_ci95"]),
        "aggregation_counts": point["counts"],
        "contrast_count": len(rows),
        "dependency_cluster_count": len({row["source_dependency_cluster_id"] for row in rows}),
        "bootstrap_valid_draws": bootstrap["valid_draws"],
        "relations_sha256": sha256(relations_path),
        "numeric_method_sha256": sha256(method_path),
        "embedding_sha256": sha256(args.embeddings),
    }
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    ledgers = hierarchy_ledgers(rows, credits)
    write_tsv(args.output_dir / "INDIVIDUAL_SCORES.tsv", individual_rows)
    write_tsv(args.output_dir / "ENDPOINT_SCORES.tsv", ledgers["endpoint"])
    write_tsv(args.output_dir / "TARGET_SUPERFAMILY_SCORES.tsv", ledgers["target_sf"])
    write_tsv(args.output_dir / "ANCHOR_SCORES.tsv", ledgers["anchor"])
    write_tsv(args.output_dir / "CLUSTER_SCORES.tsv", ledgers["cluster"])
    np.save(args.output_dir / "F_BOOTSTRAP.npy", overall_draws, allow_pickle=False)
    manifest = {}
    for path in sorted(args.output_dir.iterdir()):
        if path.is_file() and path.name != "OUTPUT_MANIFEST.json":
            manifest[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(json.dumps({"formal_interface": "F_ONLY", "files": manifest}, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
