"""Validate and score the C-ST-01 deterministic Level-0 geometry reference."""

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
import score_st01_direct as public_scorer


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
        raise ValueError(f"empty table: {path.name}")
    return rows


def require_sha(value: str, label: str) -> None:
    if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"invalid {label}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", required=True, type=Path)
    parser.add_argument("--pair-evidence", required=True, type=Path)
    parser.add_argument("--objects", required=True, type=Path)
    parser.add_argument(
        "--coordinate-root",
        type=Path,
        help="Optional immutable audit mirror containing coordinate files by locator basename.",
    )
    parser.add_argument("--numeric-method", required=True, type=Path)
    parser.add_argument("--binding", required=True, type=Path)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    parser.add_argument("--individual-scores", type=Path, help="Optional contrast-level public score table")
    args = parser.parse_args()

    numeric = json.loads(args.numeric_method.read_text(encoding="utf-8"))
    binding = json.loads(args.binding.read_text(encoding="utf-8"))

    for key, path in {"relation_table_sha256": args.relations, "objects_table_sha256": args.objects, "pair_evidence_table_sha256": args.pair_evidence, "numeric_method_sha256": args.numeric_method}.items():
        if binding[key] != sha256_file(path):
            raise ValueError(f"input binding mismatch: {key}")
    if numeric["formal_reporting_score"] != "F = unweighted macro mean of dependency-cluster F values":
        raise ValueError("unexpected ST01 F-only method")
    if float(numeric["tie_tolerance"]) != 0.0:
        raise ValueError("unexpected ST01 tie tolerance")

    objects = read_tsv(args.objects)
    if len({row["object_id"] for row in objects}) != len(objects):
        raise ValueError("duplicate object_id")
    object_map = {row["object_id"]: row for row in objects}
    for row in objects:
        require_sha(row["coordinate_sha256"], "coordinate SHA256")

    evidence = read_tsv(args.pair_evidence)
    evidence_map: dict[tuple[str, str], dict[str, str]] = {}
    symmetrization_failures = 0
    for row in evidence:
        key = tuple(sorted((row["object_left_id"], row["object_right_id"])))
        if key in evidence_map:
            raise ValueError("duplicate unordered pair evidence")
        left = float(row["tm_left_normalized"])
        right = float(row["tm_right_normalized"])
        score = float(row["structural_similarity"])
        if not all(math.isfinite(x) for x in (left, right, score)):
            raise ValueError("nonfinite structure evidence")
        if abs(score - min(left, right)) > 1e-12:
            symmetrization_failures += 1
        evidence_map[key] = row
    if symmetrization_failures:
        raise ValueError("TM-score symmetrization failure")

    relations = read_tsv(args.relations)
    if len({row["contrast_id"] for row in relations}) != len(relations):
        raise ValueError("duplicate contrast_id")
    margins: list[float] = []
    pair_scores: dict[tuple[str, str], float] = {}
    required_objects: set[str] = set()
    required_pairs: set[tuple[str, str]] = set()
    edge_directions: dict[str, set[tuple[str, str]]] = defaultdict(set)
    coordinate_files_checked = 0
    individual_rows: list[dict[str, object]] = []

    for row in relations:
        anchor = row["anchor_object_id"]
        true_target = row["true_target_object_id"]
        negative = row["negative_object_id"]
        required_objects.update((anchor, true_target, negative))
        edge_directions[row["positive_edge_id"]].add((anchor, true_target))
        positive_key = tuple(sorted((anchor, true_target)))
        negative_key = tuple(sorted((anchor, negative)))
        required_pairs.update((positive_key, negative_key))
        try:
            positive = evidence_map[positive_key]
            negative_row = evidence_map[negative_key]
        except KeyError as error:
            raise ValueError("formal relation lacks frozen pair evidence") from error
        if positive["structural_class"] != "positive" or float(positive["structural_similarity"]) < 0.45:
            raise ValueError("formal positive violates frozen geometry truth")
        if negative_row["structural_class"] != "negative" or float(negative_row["structural_similarity"]) > 0.35:
            raise ValueError("formal negative violates frozen geometry truth")
        if positive["sequence_eligible"] != "True" or negative_row["sequence_eligible"] != "True":
            raise ValueError("formal relation violates frozen sequence eligibility")
        margin = float(positive["structural_similarity"]) - float(negative_row["structural_similarity"])
        if margin < 0.10 - 1e-12:
            raise ValueError("geometry margin violates frozen separation")
        margins.append(margin)
        pair_scores[(anchor, true_target)] = float(positive["structural_similarity"])
        pair_scores[(anchor, negative)] = float(negative_row["structural_similarity"])
        tie_tolerance = float(numeric["tie_tolerance"])
        if margin > tie_tolerance:
            credit = 1.0
        elif margin < -tie_tolerance:
            credit = -1.0
        else:
            credit = 0.0
        individual_rows.append({
            "contrast_id": row["contrast_id"],
            "positive_edge_id": row["positive_edge_id"],
            "anchor_object_id": anchor,
            "true_target_object_id": true_target,
            "negative_object_id": negative,
            "source_dependency_cluster_id": row["source_dependency_cluster_id"],
            "true_structural_similarity": float(positive["structural_similarity"]),
            "negative_structural_similarity": float(negative_row["structural_similarity"]),
            "geometry_margin": margin,
            "F_credit": credit,
        })

    if not required_objects <= set(object_map):
        raise ValueError("formal relation references unknown object")
    for edge_id, directions in edge_directions.items():
        if len(directions) != 2:
            raise ValueError(f"strict positive-edge bundle is incomplete: {edge_id}")
        first, second = tuple(directions)
        if first != (second[1], second[0]):
            raise ValueError(f"positive-edge directions are not exact reversals: {edge_id}")
    for row in relations:
        anchor = row["anchor_object_id"]
        if row["source_dependency_cluster_id"] != object_map[anchor]["source_dependency_cluster_id"]:
            raise ValueError("relation dependency cluster disagrees with anchor object ledger")
    for object_id in sorted(required_objects):
        row = object_map[object_id]
        if args.coordinate_root is None:
            coordinate_path = (args.objects.parent / row["coordinate_locator"]).resolve()
        else:
            coordinate_path = (args.coordinate_root / Path(row["coordinate_locator"]).name).resolve()
        if not coordinate_path.is_file():
            raise ValueError("frozen coordinate asset is missing")
        if sha256_file(coordinate_path) != row["coordinate_sha256"]:
            raise ValueError("frozen coordinate asset hash mismatch")
        coordinate_files_checked += 1

    public_result, public_ledgers, _ = public_scorer.aggregate(
        relations,
        pair_scores,
        float(numeric["tie_tolerance"]),
        int(numeric["bootstrap"]["seed"]),
        int(numeric["bootstrap"]["draws"]),
    )
    public_by_contrast = {row["contrast_id"]: row for row in public_ledgers["contrast"]}
    if any(public_by_contrast[row["contrast_id"]]["F_credit"] != row["F_credit"] for row in individual_rows):
        raise ValueError("public scorer credit mismatch")
    f_value = public_result["F"]
    if f_value != 1.0:
        raise ValueError("definitional geometry reference did not score exactly one")

    counts = {
        "formal_objects": len(required_objects),
        "positive_edges": len(edge_directions),
        "directed_contrasts": len(relations),
        "endpoint_challenges": public_result["counts"]["endpoint_challenges"],
        "anchors": public_result["counts"]["anchors"],
        "source_dependency_clusters": public_result["counts"]["source_dependency_clusters"],
        "required_unordered_pairs": len(required_pairs),
    }
    result = {
        "constraint_id": "C-ST-01",
        "formal_interface": "F_ONLY",
        "score_transform_version": "ST01_v0.1.2_F_ONLY",
        "reference_tier": "LEVEL_0_DETERMINISTIC_RAW_INPUT_GEOMETRY",
        "reference_role": "DEFINITIONAL_WHITE_BOX_REFERENCE",
        "score_quantity": "FROZEN_MIN_BIDIRECTIONAL_LENGTH_NORMALIZED_TM_SCORE",
        "F": f_value,
        "uncertainty": {
            "status": "NOT_APPLICABLE_DEFINITIONAL_REFERENCE",
            "F_95_CI": None,
            "reason": "All formal contrasts are deterministically separated by the frozen truth thresholds; this is not an empirical comparator estimate.",
        },
        "minimum_geometry_margin": min(margins),
        "counts": counts,
        "input_identity": {
            "relations_tsv_sha256": sha256_file(args.relations),
            "pair_evidence_tsv_sha256": sha256_file(args.pair_evidence),
            "objects_tsv_sha256": sha256_file(args.objects),
            "numeric_method_sha256": sha256_file(args.numeric_method),
            "input_binding_sha256": sha256_file(args.binding),
        },
        "status": "REVIEW_RESULT",
    }
    validation = {
        "status": "PASS",
        "truth_integrity_failure_count": 0,
        "coordinate_hash_closure": {"checked": coordinate_files_checked, "failures": 0},
        "pair_symmetrization_failure_count": 0,
        "formal_relation_evidence_closure": "PASS",
        "strict_edge_bundle_closure": "PASS",
        "sequence_eligibility_closure": "PASS",
        "geometry_threshold_and_margin_closure": "PASS",
        "aggregation_hierarchy": "contrast -> endpoint -> anchor -> source dependency cluster -> unweighted cluster macro",
        "result": {"F": f_value},
        "counts": counts,
    }
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.validation.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.validation.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.individual_scores:
        args.individual_scores.parent.mkdir(parents=True, exist_ok=True)
        with args.individual_scores.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(individual_rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(individual_rows)
    print(json.dumps({"status": "PASS", "F": f_value, **counts}, sort_keys=True))


if __name__ == "__main__":
    main()
