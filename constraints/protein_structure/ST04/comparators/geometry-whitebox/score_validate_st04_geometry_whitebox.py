"""Validate and score the C-ST-04 deterministic coordinate-geometry reference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path


SCORING_DIR = Path(__file__).resolve().parents[2] / "scoring"
if str(SCORING_DIR) not in sys.path:
    sys.path.insert(0, str(SCORING_DIR))
import st04_numeric_core_v0_1_1 as numeric_core


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("empty relations table")
    return rows


def mean(values: list[float]) -> float:
    if not values:
        raise ValueError("empty aggregation unit")
    return sum(values) / len(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", required=True, type=Path)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    args = parser.parse_args()

    score = {
        "positive_field": "positive_distance", "negative_field": "negative_distance",
        "tie_tolerance": 0.0, "positive_max_A": 8.0, "negative_min_A": 12.0,
    }

    rows = read_tsv(args.relations)
    relation_ids: set[str] = set()
    queries: dict[str, dict[str, object]] = {}
    margins: list[float] = []
    wins = ties = losses = 0

    for row in rows:
        rid = row["relation_id"]
        if not rid or rid in relation_ids:
            raise ValueError("missing or duplicate relation_id")
        relation_ids.add(rid)
        if row["matching_status"] != "MATCHED_SELECTED":
            raise ValueError("nonformal matching status")
        if row["aa_exact_match"].lower() != "true" or row["target_aa_identity"] != row["negative_aa_identity"]:
            raise ValueError("exact-AA matching failure")
        if row["break_segment_eligibility"] != "SAME_FROZEN_COORDINATE_SEGMENT":
            raise ValueError("coordinate-segment closure failure")
        if row["signed_direction_class"] not in {"N_TERMINAL", "C_TERMINAL"}:
            raise ValueError("signed-direction closure failure")
        positive = float(row[score["positive_field"]])
        negative = float(row[score["negative_field"]])
        if not math.isfinite(positive) or not math.isfinite(negative):
            raise ValueError("nonfinite geometry")
        if positive > score["positive_max_A"] + 1e-12:
            raise ValueError("positive geometry threshold failure")
        if negative < score["negative_min_A"] - 1e-12:
            raise ValueError("negative geometry threshold failure")
        margin = negative - positive
        margins.append(margin)
        tau = score["tie_tolerance"]
        if margin > tau:
            credit = 1.0
            wins += 1
        elif margin < -tau:
            credit = -1.0
            losses += 1
        else:
            credit = 0.0
            ties += 1

        qid = row["query_id"]
        base = (
            row["source_dependency_cluster_id"],
            row["structure_object_id"],
            row["anchor_residue_id"],
            row["true_target_residue_id"],
        )
        q = queries.setdefault(qid, {"base": base, "rank_credit": []})
        if q["base"] != base:
            raise ValueError("query endpoint closure failure")
        q["rank_credit"].append((int(row["within_query_rank"]), credit))

    public_queries: dict[str, dict[str, str]] = {}
    public_predictions: dict[str, dict[str, float]] = {}
    query_k_counts: dict[int, int] = defaultdict(int)
    for qid, q in queries.items():
        rank_credit = sorted(q["rank_credit"])
        if [rank for rank, _ in rank_credit] != list(range(1, len(rank_credit) + 1)):
            raise ValueError(f"query rank closure failure: {qid}")
        query_credit = mean([credit for _, credit in rank_credit])
        cluster, structure, anchor, _ = q["base"]
        public_queries[qid] = {"cluster": cluster, "structure": structure, "anchor": anchor}
        public_predictions[qid] = {"query_credit": query_credit}
        query_k_counts[len(rank_credit)] += 1

    f_value, cluster_values, hierarchy = numeric_core.aggregate(
        public_queries, public_predictions, return_hierarchy=True
    )
    anchor_count = len(hierarchy["anchor"])
    structure_count = len(hierarchy["structure"])
    if f_value != 1.0 or ties or losses:
        raise ValueError("TRUTH_INTEGRITY_FAILURE: geometry reference did not score exactly one")

    counts = {
        "formal_relations": len(rows),
        "directed_queries": len(queries),
        "anchors": anchor_count,
        "formal_structures": structure_count,
        "dependency_clusters": len(cluster_values),
        "wins": wins,
        "ties": ties,
        "losses": losses,
        "query_negative_count_distribution": {str(k): query_k_counts[k] for k in sorted(query_k_counts)},
    }
    result = {
        "artifact_id": "ST04_GEOMETRY_WHITEBOX_F_ONLY_RESULT",
        "constraint_id": "C-ST-04",
        "reference_tier": "LEVEL_0_DETERMINISTIC_RAW_INPUT_GEOMETRY",
        "reference_role": "DEFINITIONAL_WHITE_BOX_REFERENCE",
        "score_quantity": "NEGATIVE_FROZEN_CA_EUCLIDEAN_DISTANCE",
        "formal_interface": "F_ONLY",
        "F": f_value,
        "uncertainty": {
            "status": "NOT_APPLICABLE_DEFINITIONAL_REFERENCE",
            "F_95_CI": None,
            "reason": "The result is guaranteed by the frozen truth separation after full integrity validation; it is not an empirical predictive-performance estimate.",
        },
        "minimum_geometry_margin_A": min(margins),
        "maximum_geometry_margin_A": max(margins),
        "counts": counts,
        "relations_sha256": sha256_file(args.relations),
        "status": "FROZEN_RESULT",
    }
    validation = {
        "artifact_id": "ST04_GEOMETRY_WHITEBOX_VALIDATION_v0.1.0_FROZEN_20261012",
        "status": "PASS",
        "truth_integrity_failure_count": 0,
        "input_closure": "PASS",
        "relation_and_query_closure": "PASS",
        "geometry_threshold_closure": "PASS",
        "exact_aa_direction_separation_segment_matching_closure": "PASS",
        "aggregation_hierarchy": "contrast -> directed query -> anchor -> structure -> dependency cluster -> unweighted cluster macro",
        "result": {"F": f_value},
        "minimum_geometry_margin_A": min(margins),
        "counts": counts,
    }
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.validation.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.validation.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "F": f_value, **counts}, sort_keys=True))


if __name__ == "__main__":
    main()
