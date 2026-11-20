#!/usr/bin/env python3
"""Public formal PS04 scorer for a complete frozen leaf ledger."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

SCORING_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCORING_DIR))
from ps04_anchor_source_unit_bootstrap import build_maps, bootstrap, point_source_unit_values
from ps04_formal_population_gate import assert_complete_population
from ps04_scoring_core import load_rows, score


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_units(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    out: dict[str, str] = {}
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            anchor, unit = row["anchor_id"], row["anchor_source_unit_bootstrap_unit"]
            previous = out.setdefault(anchor, unit)
            if previous != unit:
                raise ValueError(f"non-unique frozen uncertainty unit for {anchor}")
    if not out:
        raise ValueError("empty uncertainty membership")
    return out


def anchor_metadata(rows: list[dict], unit_by_anchor: dict[str, str]) -> tuple[list[str], list[dict], list[str]]:
    fields = ("anchor_id", "anchor_parent_id", "anchor_source_unit_id", "anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id")
    meta: dict[str, dict] = {}
    for row in rows:
        anchor = row["anchor_id"]
        value = {field: row[field] for field in fields}
        previous = meta.setdefault(anchor, value)
        if previous != value:
            raise ValueError(f"inconsistent anchor metadata for {anchor}")
    if set(meta) != set(unit_by_anchor):
        raise ValueError("leaf population and frozen uncertainty membership have different anchor sets")
    anchors = sorted(meta)
    return anchors, [meta[a] for a in anchors], [unit_by_anchor[a] for a in anchors]


def percentile(values, p: float) -> float:
    import numpy as np
    return float(np.quantile(values, p, method="linear"))


def evaluation_id(row: dict) -> str:
    identity = [row[key] for key in (
        "negative_endpoint_id", "negative_hsf_id", "negative_topology_id",
        "positive_endpoint_id", "positive_hsf_id", "anchor_id",
        "anchor_parent_id", "anchor_source_unit_id", "anchor_hsf_id",
        "anchor_topology_id", "anchor_architecture_id",
    )]
    payload = json.dumps(identity, separators=(",", ":"), ensure_ascii=False)
    return "PS04_EVAL:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_individual_scores(path: Path, rows: list[dict], unit_by_anchor: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    fields = [
        "evaluation_id", "anchor_id", "positive_endpoint_id", "negative_endpoint_id",
        "anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id",
        "anchor_source_unit_id", "bootstrap_unit_id", "positive_score",
        "negative_score", "margin", "F_credit",
    ]
    with opener(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "evaluation_id": evaluation_id(row),
                "anchor_id": row["anchor_id"],
                "positive_endpoint_id": row["positive_endpoint_id"],
                "negative_endpoint_id": row["negative_endpoint_id"],
                "anchor_hsf_id": row["anchor_hsf_id"],
                "anchor_topology_id": row["anchor_topology_id"],
                "anchor_architecture_id": row["anchor_architecture_id"],
                "anchor_source_unit_id": row["anchor_source_unit_id"],
                "bootstrap_unit_id": unit_by_anchor[row["anchor_id"]],
                "positive_score": repr(float(row["positive_score"])),
                "negative_score": repr(float(row["negative_score"])),
                "margin": repr(float(row["margin"])),
                "F_credit": repr(2.0 * float(row["value"]) - 1.0),
            })


def write_aggregation_trace(path: Path, trace: dict[str, list[dict]]) -> None:
    """Persist the compact anchor-level bridge from individual credits to F."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        # Relation-level membership is already materialized in
        # INDIVIDUAL_SCORES.tsv.gz.  Intermediate negative-side group nodes
        # would reproduce that population without adding public evidence.
        for level, nodes in trace.items():
            if level != "ANCHOR":
                continue
            for node in nodes:
                child_count = len(node["child_ids"])
                record = {
                    "node_id": node["node_id"],
                    "level": node["level"],
                    "group_keys": node["group_keys"],
                    # Leaf membership is recoverable from INDIVIDUAL_SCORES.tsv.gz
                    # using group_keys.  Repeating all child IDs at every level
                    # expands the same relation population many times.
                    "child_count": child_count,
                    "weight_rule": "equal",
                    "child_weight": (1.0 / child_count) if child_count else None,
                    "F": 2.0 * node["value"] - 1.0,
                    "population_digest": node["population_digest"],
                }
                handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")


def frozen_high_reuse(path: Path, threshold: int) -> tuple[set[tuple[str, str]], dict[str, int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    high: set[tuple[str, str]] = set()
    by_role = Counter()
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            role = {"POSITIVE_TARGET": "positive", "NEGATIVE_TARGET": "negative"}.get(row["role"])
            if role is None:
                raise ValueError(f"unsupported frozen reuse role: {row['role']}")
            if int(row["anchor_bundle_reuse_count"]) > threshold:
                high.add((role, row["endpoint_id"]))
                by_role[role] += 1
    return high, dict(sorted(by_role.items()))


def reuse_sensitivity(rows: list[dict], population_digest: str, reuse_ledger: Path, threshold: int) -> dict:
    high, high_by_role = frozen_high_reuse(reuse_ledger, threshold)
    retained = [row for row in rows if ("positive", row["positive_endpoint_id"]) not in high and ("negative", row["negative_endpoint_id"]) not in high]
    result = {"kind": "LEAVE_HIGH_REUSE_ENDPOINT_OUT", "threshold": threshold,
              "reuse_definition": "frozen_anchor_bundle_reuse_count", "endpoint_reuse_ledger_sha256": sha(reuse_ledger),
              "high_reuse_endpoint_count": len(high), "high_reuse_endpoint_count_by_role": high_by_role, "input_leaf_count": len(rows),
              "retained_leaf_count": len(retained), "bootstrap": "NOT_RUN_SENSITIVITY_POINT_ESTIMATE_ONLY"}
    if not retained:
        return result | {"status": "NO_RETAINED_LEAVES"}
    _, raw = score(retained, population_digest + ":reuse_gt_" + str(threshold))
    return result | {"status": "COMPLETE", "F": 2.0 * raw - 1.0}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binding", type=Path, required=True)
    ap.add_argument("--formal-relations", type=Path, required=True)
    ap.add_argument("--leaf-ledger", type=Path, required=True)
    ap.add_argument("--uncertainty-membership", type=Path, required=True)
    ap.add_argument("--endpoint-reuse-ledger", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--individual-scores", type=Path)
    args = ap.parse_args()
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    refset_paths = {"formal_relation_ledger_sha256": args.formal_relations,
                    "uncertainty_membership_sha256": args.uncertainty_membership,
                    "endpoint_reuse_ledger_sha256": args.endpoint_reuse_ledger}
    refset_identity = {
        key: {"actual_sha256": sha(path), "registered_sha256": binding[key],
              "registered_byte_identity": sha(path) == binding[key]}
        for key, path in refset_paths.items()
    }
    registered_refset_byte_identity = all(
        value["registered_byte_identity"] for value in refset_identity.values()
    )
    population = assert_complete_population(args.formal_relations, args.leaf_ledger)
    rows = load_rows(args.leaf_ledger, float(binding["delta_num"]))
    trace, raw = score(rows, population["leaf_ledger_sha256"])
    units = read_units(args.uncertainty_membership)
    individual_scores = args.individual_scores or args.output_dir / "INDIVIDUAL_SCORES.tsv.gz"
    write_individual_scores(individual_scores, rows, units)
    aggregation_trace = args.output_dir / "AGGREGATION_TRACE.jsonl.gz"
    write_aggregation_trace(aggregation_trace, trace)
    anchors, meta, unit_labels = anchor_metadata(rows, units)
    maps = build_maps(meta, unit_labels)
    values = {node["group_keys"]["anchor_id"]: node["value"] for node in trace["ANCHOR"]}
    anchor_values = [values[a] for a in anchors]
    source_unit_values = point_source_unit_values(anchor_values, maps)
    # With every frozen source unit represented exactly once, this must reconstruct
    # the complete HSF -> topology -> architecture point hierarchy.
    import numpy as np
    from ps04_anchor_source_unit_bootstrap import vectorized_draws
    point_reconstruction = float(vectorized_draws(source_unit_values, maps, np.ones((1, maps["n_bootstrap"])))[0])
    if abs(point_reconstruction - raw) > 1e-12:
        raise ValueError(f"source-unit point reconstruction drift: {point_reconstruction} vs {raw}")
    primary = binding["primary_bootstrap"]
    draws = bootstrap(source_unit_values, maps, int(primary["draws"]), int(primary["seed"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sensitivity = reuse_sensitivity(rows, population["leaf_ledger_sha256"], args.endpoint_reuse_ledger, int(binding["reuse_sensitivity"]["threshold"]))
    result = {"artifact_id": "PS04_FORMAL_RESULT_v0.1.2", "status": "FORMAL_EXECUTION_COMPLETE",
              "refset_release_id": binding["refset_release_id"], "binding_sha256": sha(args.binding),
              "refset_input_identity": refset_identity,
              "registered_refset_byte_identity": registered_refset_byte_identity,
              "input_refset_status": "REGISTERED_RELEASE_BYTES" if registered_refset_byte_identity else "SEMANTICALLY_VALID_COMPATIBLE_INPUT",
              "population_gate": population, "delta_num": binding["delta_num"],
              "bootstrap_unit_count": maps["n_bootstrap"], "draws": int(primary["draws"]), "rng": primary["rng"], "seed": int(primary["seed"]),
              "formal_interface": "F_ONLY", "score_transform_version": "PS04_v1.0.2_F_ONLY",
              "F": 2.0 * raw - 1.0,
              "F_95_CI": [2.0 * percentile(draws, .025) - 1.0, 2.0 * percentile(draws, .975) - 1.0],
              "point_reconstruction_F": 2.0 * point_reconstruction - 1.0,
              "reuse_sensitivity": sensitivity,
              "coverage": {"formal_evaluation_units": len(rows), "scored_evaluation_units": len(rows), "status": "COMPLETE"},
              "individual_scores": {"path": str(individual_scores), "sha256": sha(individual_scores)},
              "aggregation_trace": {"path": str(aggregation_trace), "sha256": sha(aggregation_trace),
                                    "leaf_level": "INDIVIDUAL_SCORES.tsv.gz", "included_levels": ["ANCHOR"]},
              "code_sha256": {
                  "score_ps04_formal.py": sha(Path(__file__)),
                  "ps04_scoring_core.py": sha(SCORING_DIR / "ps04_scoring_core.py"),
                  "ps04_anchor_source_unit_bootstrap.py": sha(SCORING_DIR / "ps04_anchor_source_unit_bootstrap.py"),
                  "ps04_formal_population_gate.py": sha(Path(__file__).with_name("ps04_formal_population_gate.py")),
              }}
    (args.output_dir / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()


