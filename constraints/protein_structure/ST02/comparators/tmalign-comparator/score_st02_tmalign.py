from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scoring"))
from st02_numeric_core import aggregate_hierarchy, bayesian_cluster_bootstrap
from score_st02 import hierarchy_ledgers, stringify, write_tsv


SCORER_ID = "ST02_TMALIGN_COMPARATOR_SCORER_v0.1.3_F_ONLY"
HERE = Path(__file__).resolve().parent
CONSTRAINT = HERE.parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-ledger", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--individual-output", type=Path)
    args = parser.parse_args()

    ledger = load_json(args.pair_ledger)

    relations_path = CONSTRAINT / "refset" / "relations.tsv"
    numeric_path = CONSTRAINT / "scoring" / "ST02_NUMERIC_METHOD_v0.1.3_F_ONLY.json"
    if not relations_path.is_file() or not numeric_path.is_file():
        raise ValueError("public RefSet or Numeric Method is unavailable")
    numeric = load_json(numeric_path)
    if numeric.get("status") != "FROZEN_NUMERIC_METHOD":
        raise ValueError("Numeric Method not frozen")
    with relations_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    pair_scores = {}
    pair_rows = ledger.get("pairs", [])
    for row in pair_rows:
        key = tuple(sorted((row["object_a"], row["object_b"])))
        score = float(row["symmetric_tm_score"])
        expected = (float(row["tm_score_chain_1_normalized"]) + float(row["tm_score_chain_2_normalized"])) / 2.0
        if key in pair_scores or not math.isfinite(score) or score != expected:
            raise ValueError("duplicate, nonfinite, or inconsistent pair score")
        pair_scores[key] = score
    needed = set()
    for row in rows:
        needed.add(tuple(sorted((row["anchor_object_id"], row["true_target_object_id"]))))
        needed.add(tuple(sorted((row["anchor_object_id"], row["negative_object_id"]))))
    if set(pair_scores) != needed or ledger.get("pair_count") != len(pair_rows) or ledger.get("required_pair_count") != len(needed):
        raise ValueError("pair coverage closure failure")

    credits = []
    individual = []
    for row in rows:
        positive = pair_scores[tuple(sorted((row["anchor_object_id"], row["true_target_object_id"])))]
        negative = pair_scores[tuple(sorted((row["anchor_object_id"], row["negative_object_id"])))]
        delta = positive - negative
        credit = 1.0 if delta > 0 else -1.0 if delta < 0 else 0.0
        credits.append(credit)
        individual.append({
            "contrast_id": row["contrast_id"],
            "direction_id": row["direction_id"],
            "anchor_object_id": row["anchor_object_id"],
            "true_target_object_id": row["true_target_object_id"],
            "negative_object_id": row["negative_object_id"],
            "source_dependency_cluster_id": row["source_dependency_cluster_id"],
            "topology_id": row["topology_id"],
            "positive_tm_score": positive,
            "negative_tm_score": negative,
            "raw_margin": delta,
            "F_credit": credit,
        })
    point = aggregate_hierarchy(rows, credits)
    uncertainty = bayesian_cluster_bootstrap(
        rows, credits, draws=int(numeric["uncertainty"]["draws"]), seed=int(numeric["uncertainty"]["seed"]),
        return_overall_draws=True,
    )
    bootstrap_draws = uncertainty.pop("overall_draws")
    result = {
        "artifact_id": "ST02_TMALIGN_COMPARATOR_F_ONLY_RESULT", "status": "FORMAL_COMPARATOR_RESULT",
        "formal_interface": "F_ONLY", "score_transform_version": "ST02_v0.1.3_F_ONLY",
        "scorer_id": SCORER_ID, "F": point["F"], "F_95_CI": uncertainty["overall_ci95"],
        "source_sf_topology_F": stringify(point["source_sf_topology"]),
        "source_sf_topology_F_95_CI": stringify(uncertainty["source_sf_topology_ci95"]),
        "topology_F": stringify(point["topology"]),
        "topology_F_95_CI": stringify(uncertainty["topology_ci95"]), "pair_count": len(pair_scores),
        "contrast_count": len(rows), "aggregation_counts": point["counts"],
        "pair_ledger_sha256": sha256(args.pair_ledger),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    individual_output = args.individual_output or args.output.with_name("INDIVIDUAL_SCORES.tsv")
    individual_output.parent.mkdir(parents=True, exist_ok=True)
    with individual_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(individual[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(individual)
    ledgers = hierarchy_ledgers(rows, credits)
    for filename, key in (
        ("ENDPOINT_SCORES.tsv", "endpoint"),
        ("TARGET_SUPERFAMILY_SCORES.tsv", "target_sf"),
        ("ANCHOR_SCORES.tsv", "anchor"),
        ("CLUSTER_SCORES.tsv", "cluster"),
    ):
        write_tsv(args.output.parent / filename, ledgers[key])
    np.save(args.output.parent / "F_BOOTSTRAP.npy", bootstrap_draws, allow_pickle=False)
    outputs = [
        args.output,
        individual_output,
        *(args.output.parent / name for name in (
            "ENDPOINT_SCORES.tsv", "TARGET_SUPERFAMILY_SCORES.tsv",
            "ANCHOR_SCORES.tsv", "CLUSTER_SCORES.tsv", "F_BOOTSTRAP.npy",
        )),
    ]
    manifest = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in outputs
    }
    (args.output.parent / "RESULT_MANIFEST.json").write_text(
        json.dumps({"formal_interface": "F_ONLY", "files": manifest}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
