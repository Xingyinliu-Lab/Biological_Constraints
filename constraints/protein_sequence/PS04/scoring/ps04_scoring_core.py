#!/usr/bin/env python3
"""Schema-bound PS04 aggregation kernel.

Consumes a frozen, complete leaf-score TSV.  It does not construct a RefSet,
select endpoints, run a representation model, or change the raw estimand.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

LEVELS = [
    "NEGATIVE_ENDPOINT", "NEGATIVE_HSF", "NEGATIVE_TOPOLOGY",
    "POSITIVE_ENDPOINT", "POSITIVE_HSF", "ANCHOR", "SOURCE_PARENT",
    "SOURCE_UNIT_WITHIN_ANCHOR_HSF_TOPOLOGY", "ANCHOR_HSF",
    "ANCHOR_TOPOLOGY", "ANCHOR_ARCHITECTURE", "OVERALL",
]
REQUIRED = [
    "negative_endpoint_id", "negative_hsf_id", "negative_topology_id",
    "positive_endpoint_id", "positive_hsf_id", "anchor_id",
    "anchor_parent_id", "anchor_source_unit_id", "anchor_hsf_id",
    "anchor_topology_id", "anchor_architecture_id", "positive_score",
    "negative_score",
]

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)

def node_id(level: str, keys: list[str]) -> str:
    return f"{level}:{hashlib.sha256(canonical_json([level, keys]).encode()).hexdigest()}"

def validate_contract(root: Path) -> dict:
    metric = json.loads((root / "contracts" / "MetricRegistry.json").read_text(encoding="utf-8-sig"))
    graph = json.loads((root / "contracts" / "AggregationGraph.json").read_text(encoding="utf-8-sig"))
    result = json.loads((root / "contracts" / "ResultContract.json").read_text(encoding="utf-8-sig"))
    entries = {x["metric_id"]: x for x in metric["contract"]["entries"]}
    formal = entries.get("F_PS04", {})
    if formal.get("role") != "PRIMARY" or formal.get("range") != [-1, 1]:
        raise ValueError("MetricRegistry must declare F_PS04 as the primary formal metric")
    if graph["contract"]["levels"] != LEVELS or result["contract"]["primary_metrics"] != ["F_PS04"]:
        raise ValueError("Aggregation graph or primary metric declaration drift")
    return {
        "metric_registry_sha256": digest(root / "contracts" / "MetricRegistry.json"),
        "aggregation_graph_sha256": digest(root / "contracts" / "AggregationGraph.json"),
        "result_contract_sha256": digest(root / "contracts" / "ResultContract.json"),
    }

def load_rows(path: Path, delta: float) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(REQUIRED) - set(reader.fieldnames):
            missing = sorted(set(REQUIRED) - set(reader.fieldnames or []))
            raise ValueError(f"leaf ledger lacks required fields: {missing}")
        rows = []
        seen = set()
        for idx, row in enumerate(reader, start=2):
            identity = tuple(row[k] for k in REQUIRED[:-2])
            if identity in seen:
                raise ValueError(f"duplicate leaf identity at input line {idx}")
            seen.add(identity)
            try:
                positive, negative = float(row["positive_score"]), float(row["negative_score"])
            except ValueError as exc:
                raise ValueError(f"non-numeric score at input line {idx}") from exc
            if not math.isfinite(positive) or not math.isfinite(negative):
                raise ValueError(f"non-finite score at input line {idx}")
            margin = positive - negative
            credit = 1.0 if margin > delta else 0.0 if margin < -delta else 0.5
            row["margin"] = margin
            row["value"] = credit
            rows.append(row)
    if not rows:
        raise ValueError("leaf ledger is empty")
    return rows

def assert_constant(rows: list[dict], key: str, attributes: list[str]) -> None:
    values = {}
    for row in rows:
        signature = tuple(row[a] for a in attributes)
        if row[key] in values and values[row[key]] != signature:
            raise ValueError(f"inconsistent identity mapping for {key}={row[key]}")
        values[row[key]] = signature

def aggregate(children: list[dict], level: str, fields: list[str], population_digest: str) -> list[dict]:
    groups = defaultdict(list)
    for child in children:
        groups[tuple(child["meta"][f] for f in fields)].append(child)
    nodes = []
    for key in sorted(groups):
        group = groups[key]
        ids = [x["node_id"] for x in group]
        values = [x["value"] for x in group]
        nodes.append({
            "node_id": node_id(level, list(key)), "level": level,
            "group_keys": dict(zip(fields, key)), "child_ids": ids,
            "child_values": values, "weights": [1.0 / len(group)] * len(group),
            "value": math.fsum(values) / len(values), "population_digest": population_digest,
            "meta": group[0]["meta"],
        })
    return nodes

def score(rows: list[dict], population_digest: str) -> tuple[dict[str, list[dict]], float]:
    for key, attrs in [
        ("anchor_id", ["anchor_parent_id", "anchor_source_unit_id", "anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id"]),
        ("positive_endpoint_id", ["positive_hsf_id"]),
        ("negative_endpoint_id", ["negative_hsf_id", "negative_topology_id"]),
    ]:
        assert_constant(rows, key, attrs)
    trace = {"NEGATIVE_ENDPOINT": []}
    for row in rows:
        meta = {k: row[k] for k in REQUIRED[:-2]}
        trace["NEGATIVE_ENDPOINT"].append({
            "node_id": node_id("NEGATIVE_ENDPOINT", [row["negative_endpoint_id"], row["positive_endpoint_id"], row["anchor_id"]]),
            "level": "NEGATIVE_ENDPOINT", "group_keys": {"negative_endpoint_id": row["negative_endpoint_id"], "positive_endpoint_id": row["positive_endpoint_id"], "anchor_id": row["anchor_id"]},
            "child_ids": [], "child_values": [], "weights": [], "value": row["value"],
            "population_digest": population_digest, "meta": meta,
        })
    transitions = [
        ("NEGATIVE_HSF", ["anchor_id", "positive_endpoint_id", "negative_topology_id", "negative_hsf_id"]),
        ("NEGATIVE_TOPOLOGY", ["anchor_id", "positive_endpoint_id", "negative_topology_id"]),
        ("POSITIVE_ENDPOINT", ["anchor_id", "positive_endpoint_id"]),
        ("POSITIVE_HSF", ["anchor_id", "positive_hsf_id"]),
        ("ANCHOR", ["anchor_id"]),
        ("SOURCE_PARENT", ["anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id", "anchor_parent_id"]),
        ("SOURCE_UNIT_WITHIN_ANCHOR_HSF_TOPOLOGY", ["anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id", "anchor_source_unit_id"]),
        ("ANCHOR_HSF", ["anchor_architecture_id", "anchor_topology_id", "anchor_hsf_id"]),
        ("ANCHOR_TOPOLOGY", ["anchor_architecture_id", "anchor_topology_id"]),
        ("ANCHOR_ARCHITECTURE", ["anchor_architecture_id"]),
        ("OVERALL", []),
    ]
    children = trace["NEGATIVE_ENDPOINT"]
    for level, fields in transitions:
        children = aggregate(children, level, fields, population_digest)
        trace[level] = children
    if len(trace["OVERALL"]) != 1:
        raise ValueError("PS04 graph must yield exactly one OVERALL node")
    return trace, trace["OVERALL"][0]["value"]

def main() -> None:
    ap = argparse.ArgumentParser(description="Score a complete frozen PS04 leaf ledger.")
    ap.add_argument("--contract-root", type=Path, required=True)
    ap.add_argument("--leaf-ledger", type=Path, required=True)
    ap.add_argument("--delta-num", type=float, required=True, help="Frozen numeric preflight tolerance; supplied explicitly.")
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    if args.delta_num < 0 or not math.isfinite(args.delta_num):
        raise ValueError("--delta-num must be finite and non-negative")
    binding = validate_contract(args.contract_root)
    rows = load_rows(args.leaf_ledger, args.delta_num)
    population = digest(args.leaf_ledger)
    trace, raw = score(rows, population)
    signed = 2.0 * raw - 1.0
    args.output_dir.mkdir(parents=True, exist_ok=True)
    flat = []
    for level in LEVELS:
        for node in trace[level]:
            public = dict(node, meta=None)
            public["F"] = 2.0 * float(public.pop("value")) - 1.0
            public["child_F"] = [2.0 * float(x) - 1.0 for x in public.pop("child_values")]
            flat.append(public)
    (args.output_dir / "AGGREGATION_TRACE.jsonl").write_text("\n".join(canonical_json(x) for x in flat) + "\n", encoding="utf-8")
    summary = {
        "kind": "PS04_SCORE_KERNEL_RESULT", "formal_result": False,
        "input_leaf_ledger_sha256": population, "delta_num": args.delta_num,
        "leaf_count": len(rows), "level_counts": {level: len(trace[level]) for level in LEVELS},
        "F": {"metric_id": "F_PS04", "value": signed, "range": [-1, 1], "neutral": 0},
        "contract_binding": binding,
    }
    (args.output_dir / "SCORE_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

if __name__ == "__main__":
    main()
