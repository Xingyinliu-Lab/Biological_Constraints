"""Export ID-level measurements from a completed frozen ST04 geometry reference run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", required=True, type=Path)
    parser.add_argument("--formal-result", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    result = json.loads(args.formal_result.read_text(encoding="utf-8"))
    if result.get("status") not in {"PASS", "FROZEN_RESULT"} or result.get("F") != 1.0:
        raise ValueError("the frozen geometry comparator must pass before exporting individual scores")

    with args.relations.open(encoding="utf-8", newline="") as handle:
        relations = list(csv.DictReader(handle, delimiter="\t"))
    if not relations:
        raise ValueError("empty relations table")

    output_rows = []
    relation_ids = set()
    for row in relations:
        relation_id = row["relation_id"]
        if not relation_id or relation_id in relation_ids:
            raise ValueError("missing or duplicate relation_id")
        relation_ids.add(relation_id)
        positive = float(row["positive_distance"])
        negative = float(row["negative_distance"])
        if not math.isfinite(positive) or not math.isfinite(negative):
            raise ValueError("nonfinite geometry")
        margin = negative - positive
        credit = 1.0 if margin > 0.0 else (-1.0 if margin < 0.0 else 0.0)
        output_rows.append(
            {
                "relation_id": relation_id,
                "query_id": row["query_id"],
                "source_dependency_cluster_id": row["source_dependency_cluster_id"],
                "structure_object_id": row["structure_object_id"],
                "anchor_residue_id": row["anchor_residue_id"],
                "true_target_residue_id": row["true_target_residue_id"],
                "negative_residue_id": row["negative_residue_id"],
                "positive_distance_A": positive,
                "negative_distance_A": negative,
                "geometry_margin_A": margin,
                "F_credit": credit,
            }
        )
    expected_relations = result.get("counts", {}).get("formal_relations", len(output_rows))
    if len(output_rows) != expected_relations:
        raise ValueError("individual-score row count does not match the frozen result")
    if any(row["F_credit"] != 1.0 for row in output_rows):
        raise ValueError("geometry truth-integrity failure")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
