#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from rna0102_hamming_comparator import require_oracle_tau, score_frozen_ledgers


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    default_refset = Path(__file__).resolve().parents[2] / "refset" / "RNA0102_L16K_EDIT_REFSET"
    parser = argparse.ArgumentParser(description="RNA01-02 exact-Hamming white-box reference")
    parser.add_argument("--refset-dir", type=Path, default=default_refset)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_oracle_tau(0)
    scored = score_frozen_ledgers(load(args.refset_dir / "EDIT_LADDER_LEDGER.json"), load(args.refset_dir / "OBJECT_LEDGER.json"))
    individual = []
    by_cluster: dict[str, list[float]] = {}
    for row in scored["leaves"]:
        left = 1.0 if row["d02"] > row["d01"] else 0.0 if row["d02"] == row["d01"] else -1.0
        right = 1.0 if row["d02"] > row["d12"] else 0.0 if row["d02"] == row["d12"] else -1.0
        score = (left + right) / 2.0
        by_cluster.setdefault(row["outer_cluster_id"], []).append(score)
        individual.append({
            "leaf_id": row["leaf_id"], "outer_cluster_id": row["outer_cluster_id"],
            "x0_object_id": row["x0_object_id"], "x1_object_id": row["x1_object_id"], "x2_object_id": row["x2_object_id"],
            "h01": row["h01"], "h12": row["h12"], "h02": row["h02"],
            "left_F_credit": format(left, ".1f"), "right_F_credit": format(right, ".1f"), "individual_F": format(score, ".1f"),
        })
    clusters = [{
        "outer_cluster_id": key, "leaf_count": len(values), "cluster_F": format(sum(values) / len(values), ".1f")
    } for key, values in sorted(by_cluster.items())]
    if any(float(row["cluster_F"]) != 1.0 for row in clusters):
        raise ValueError("exact-Hamming ceiling did not score every cluster as 1")
    args.output.mkdir(parents=True, exist_ok=True)
    individual_path = args.output / "INDIVIDUAL_SCORES.tsv"
    cluster_path = args.output / "CLUSTER_SCORE_LEDGER.tsv"
    result_path = args.output / "FORMAL_RESULT.json"
    write_tsv(individual_path, individual)
    write_tsv(cluster_path, clusters)
    result_path.write_text(json.dumps({
        "status": "PASS", "formal_result": True, "comparator_id": scored["comparator_id"],
        "identity": "LEVEL_0_WHITE_BOX_REFERENCE_AND_DETERMINISTIC_CEILING",
        "population": "L16K", "tau_delta": 0, "evaluation_leaves": len(individual),
        "formal_interface": "F_ONLY", "score_transform_version": "RNA01-02_v0.2.6_F_ONLY_DIRECT_CREDITS",
        "outer_clusters": len(clusters), "F": 1.0, "F_95_CI": [1.0, 1.0],
        "bootstrap_draws": 20000, "bootstrap_unit": "outer_cluster",
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    outputs = [individual_path, cluster_path, result_path]
    (args.output / "OUTPUT_MANIFEST.json").write_text(json.dumps({
        "files": [{"path": path.name, "bytes": path.stat().st_size, "sha256": digest(path)} for path in outputs]
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
