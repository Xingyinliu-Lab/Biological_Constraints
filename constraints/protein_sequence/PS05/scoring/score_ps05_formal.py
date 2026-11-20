#!/usr/bin/env python3
"""Schema-bound and population-closed PS05 aggregation kernel."""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math
from collections import defaultdict
from pathlib import Path

LEVELS = ["NEGATIVE_ENDPOINT", "NEGATIVE_HSF", "POSITIVE_ENDPOINT", "POSITIVE_HSF", "ANCHOR", "SOURCE_PARENT", "ANCHOR_HSF", "ANCHOR_TOPOLOGY", "ANCHOR_ARCHITECTURE", "OVERALL"]
IDENTITY_FIELDS = ["negative_endpoint_id", "negative_hsf_id", "positive_endpoint_id", "positive_hsf_id", "anchor_id", "anchor_parent_id", "anchor_source_unit_id", "anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id"]
PROVENANCE_FIELDS = ["assay_id", "planned_bundle_ref", "interface_coordinate_ref", "execution_attempt_ref", "measurement_contract_ref"]
SCORE_FIELDS = ["positive_score", "negative_score"]
REQUIRED = IDENTITY_FIELDS + PROVENANCE_FIELDS + SCORE_FIELDS

def digest(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def canonical_json(value: object) -> str: return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
def canonical_bytes(value: object) -> bytes: return canonical_json(value).encode("utf-8")
def node_id(level: str, keys: list[str]) -> str: return f"{level}:{hashlib.sha256(canonical_json([level, keys]).encode()).hexdigest()}"
def load_json(path: Path) -> dict: return json.loads(path.read_text(encoding="utf-8-sig"))
def identity_records(rows: list[dict]) -> list[list[str]]: return sorted([[r[k] for k in IDENTITY_FIELDS] for r in rows])
def identity_digest(records: list[list[str]]) -> str: return hashlib.sha256(canonical_bytes(records)).hexdigest()

def validate_contract(root: Path, binding_path: Path) -> dict:
    binding = load_json(binding_path)
    if binding.get("constraint_id") != "C-PS-05": raise ValueError("contract binding constraint_id drift")
    contract_dir = root / "contracts" if (root / "contracts").is_dir() else root
    manifest_path = root / "MANIFEST.json"
    actual_manifest = digest(manifest_path) if manifest_path.is_file() else None
    expected_manifest = binding.get("contract_manifest_sha256")
    if expected_manifest is not None and expected_manifest != actual_manifest:
        raise ValueError("exact contract manifest binding failed")
    refs = binding.get("contract_refs", {})
    for artifact in ("MetricRegistry", "AggregationGraph", "ResultContract"):
        if refs.get(artifact, {}).get("sha256") != digest(contract_dir / f"{artifact}.json"):
            raise ValueError(f"exact contract binding failed for {artifact}")
    relation_path = root / "FORMAL_RELATION_LEDGER.tsv"
    expected_relation = binding.get("refset_relation_ledger_sha256")
    actual_relation = digest(relation_path)
    metric, graph, result = (load_json(contract_dir / f"{x}.json") for x in ("MetricRegistry", "AggregationGraph", "ResultContract"))
    entries = {x["metric_id"]: x for x in metric["contract"]["entries"]}
    signed = entries.get("F_PS05", {})
    if signed.get("role") not in {"PRIMARY", "DISPLAY_ONLY"}: raise ValueError("F_PS05 metric role drift")
    if graph["contract"]["levels"] != LEVELS: raise ValueError("aggregation graph drift")
    return {"contract_binding_sha256": digest(binding_path), "contract_manifest_sha256": actual_manifest, "metric_registry_sha256": digest(contract_dir/"MetricRegistry.json"), "aggregation_graph_sha256": digest(contract_dir/"AggregationGraph.json"), "result_contract_sha256": digest(contract_dir/"ResultContract.json"), "actual_relation_ledger_sha256": actual_relation, "registered_relation_ledger_sha256": expected_relation, "registered_refset_bytes_match": expected_relation == actual_relation if expected_relation is not None else None}

def load_rows(path: Path, delta: float) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(REQUIRED) - set(reader.fieldnames or []))
        if missing: raise ValueError(f"leaf ledger lacks required fields: {missing}")
        rows, seen = [], set()
        for idx, row in enumerate(reader, start=2):
            identity = tuple(row[k] for k in IDENTITY_FIELDS)
            if identity in seen: raise ValueError(f"duplicate leaf identity at input line {idx}")
            seen.add(identity)
            if any(not row[k] for k in PROVENANCE_FIELDS): raise ValueError(f"missing required provenance at input line {idx}")
            try: positive, negative = float(row["positive_score"]), float(row["negative_score"])
            except ValueError as exc: raise ValueError(f"non-numeric score at input line {idx}") from exc
            if not math.isfinite(positive) or not math.isfinite(negative): raise ValueError(f"non-finite score at input line {idx}")
            margin = positive - negative
            row.update(margin=margin, value=1.0 if margin > delta else -1.0 if margin < -delta else 0.0)
            rows.append(row)
    if not rows: raise ValueError("leaf ledger is empty")
    return rows

def validate_planned_population(path: Path, rows: list[dict]) -> dict:
    planned = load_json(path)
    if planned.get("kind") != "PS05_PLANNED_SCOREABLE_POPULATION" or planned.get("constraint_id") != "C-PS-05": raise ValueError("planned population identity drift")
    if planned.get("identity_fields") != IDENTITY_FIELDS: raise ValueError("planned population identity field drift")
    expected = planned.get("leaf_identities")
    if not isinstance(expected, list) or any(not isinstance(x, list) or len(x) != len(IDENTITY_FIELDS) for x in expected): raise ValueError("invalid planned population identities")
    expected = sorted(expected)
    if len(expected) != len({tuple(x) for x in expected}): raise ValueError("duplicate planned population identity")
    frozen_identity_hash = identity_digest(expected)
    if planned.get("scoreable_population_identity_sha256") != frozen_identity_hash: raise ValueError("planned population identity digest drift")
    if identity_records(rows) != expected: raise ValueError("observed leaf identities do not exactly equal frozen planned population")
    if {r["assay_id"] for r in rows} != {planned.get("assay_id")}: raise ValueError("assay_id does not match planned population")
    if any(r["measurement_contract_ref"] != planned.get("measurement_contract_ref") for r in rows): raise ValueError("measurement_contract_ref does not match planned population")
    return {"planned_scoreable_population_sha256": digest(path), "scoreable_population_identity_sha256": frozen_identity_hash, "planned_leaf_count": len(expected)}

def validate_relation_population(root: Path, planned_path: Path) -> dict:
    relation_path = root / "FORMAL_RELATION_LEDGER.tsv"
    planned = load_json(planned_path)
    fields = planned.get("identity_fields")
    identities = planned.get("leaf_identities")
    if fields != IDENTITY_FIELDS or not isinstance(identities, list):
        raise ValueError("invalid planned population for relation validation")
    planned_rows = [dict(zip(fields, values)) for values in identities]
    planned_relations = {
        (
            row["anchor_id"], row["positive_endpoint_id"], row["negative_endpoint_id"],
            row["anchor_hsf_id"], row["positive_hsf_id"], row["negative_hsf_id"],
            row["anchor_topology_id"], row["anchor_parent_id"],
        )
        for row in planned_rows
    }
    with relation_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "anchor_domain_id", "positive_domain_id", "negative_domain_id",
            "anchor_hsf_id", "positive_hsf_id", "negative_hsf_id",
            "topology_id", "anchor_parent_protein_id",
        }
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"relation ledger lacks required fields: {missing}")
        observed_rows = [
            (
                row["anchor_domain_id"], row["positive_domain_id"], row["negative_domain_id"],
                row["anchor_hsf_id"], row["positive_hsf_id"], row["negative_hsf_id"],
                row["topology_id"], row["anchor_parent_protein_id"],
            )
            for row in reader
        ]
    observed = set(observed_rows)
    if len(observed_rows) != len(observed) or len(observed) != len(planned_relations) or observed != planned_relations:
        raise ValueError("relation ledger is not semantically identical to the planned population")
    return {"relation_semantic_validation": "PASS", "relation_count": len(observed)}

def assert_constant(rows: list[dict], key: str, attributes: list[str]) -> None:
    values = {}
    for row in rows:
        signature = tuple(row[a] for a in attributes)
        if row[key] in values and values[row[key]] != signature: raise ValueError(f"inconsistent identity mapping for {key}={row[key]}")
        values[row[key]] = signature

def assert_parent_to_source_unit(rows: list[dict]) -> None:
    values = {}
    for row in rows:
        key = tuple(row[k] for k in ("anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id", "anchor_parent_id"))
        unit = row["anchor_source_unit_id"]
        if key in values and values[key] != unit: raise ValueError("same parent maps to multiple source units within biological stratum")
        values[key] = unit

def aggregate(children: list[dict], level: str, fields: list[str], population_identity: str) -> list[dict]:
    groups = defaultdict(list)
    for child in children: groups[tuple(child["meta"][f] for f in fields)].append(child)
    nodes = []
    for key in sorted(groups):
        group = groups[key]; values = [x["value"] for x in group]
        nodes.append({"node_id": node_id(level, list(key)), "level": level, "group_keys": dict(zip(fields, key)), "child_ids": [x["node_id"] for x in group], "child_values": values, "weights": [1/len(group)]*len(group), "value": math.fsum(values)/len(values), "scoreable_population_identity_sha256": population_identity, "meta": group[0]["meta"]})
    return nodes

def score(rows: list[dict], population_identity: str) -> tuple[dict[str, list[dict]], float]:
    for key, attrs in [("anchor_id", ["anchor_parent_id", "anchor_source_unit_id", "anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id"]), ("positive_endpoint_id", ["positive_hsf_id"]), ("negative_endpoint_id", ["negative_hsf_id"])]: assert_constant(rows, key, attrs)
    assert_parent_to_source_unit(rows)
    trace = {"NEGATIVE_ENDPOINT": []}
    for row in rows:
        meta = {k: row[k] for k in IDENTITY_FIELDS}
        trace["NEGATIVE_ENDPOINT"].append({"node_id": node_id("NEGATIVE_ENDPOINT", [row["negative_endpoint_id"], row["positive_endpoint_id"], row["anchor_id"]]), "level": "NEGATIVE_ENDPOINT", "group_keys": {"negative_endpoint_id":row["negative_endpoint_id"],"positive_endpoint_id":row["positive_endpoint_id"],"anchor_id":row["anchor_id"]}, "child_ids":[], "child_values":[], "weights":[], "value":row["value"], "scoreable_population_identity_sha256":population_identity, "meta":meta})
    transitions = [("NEGATIVE_HSF", ["anchor_id", "positive_endpoint_id", "negative_hsf_id"]), ("POSITIVE_ENDPOINT", ["anchor_id", "positive_endpoint_id"]), ("POSITIVE_HSF", ["anchor_id", "positive_hsf_id"]), ("ANCHOR", ["anchor_id"]), ("SOURCE_PARENT", ["anchor_hsf_id", "anchor_topology_id", "anchor_architecture_id", "anchor_parent_id"]), ("ANCHOR_HSF", ["anchor_architecture_id", "anchor_topology_id", "anchor_hsf_id"]), ("ANCHOR_TOPOLOGY", ["anchor_architecture_id", "anchor_topology_id"]), ("ANCHOR_ARCHITECTURE", ["anchor_architecture_id"]), ("OVERALL", [])]
    children = trace["NEGATIVE_ENDPOINT"]
    for level, fields in transitions: children = aggregate(children, level, fields, population_identity); trace[level] = children
    if len(trace["OVERALL"]) != 1: raise ValueError("PS05 graph must yield exactly one OVERALL node")
    return trace, trace["OVERALL"][0]["value"]

def write_individual(path: Path, rows: list[dict]) -> None:
    fields = ["evaluation_id", *IDENTITY_FIELDS, "positive_score", "negative_score", "margin", "F_credit"]
    opener = gzip.open if path.suffix == ".gz" else open
    path.parent.mkdir(parents=True, exist_ok=True)
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            identity = [row[key] for key in IDENTITY_FIELDS]
            writer.writerow({
                "evaluation_id": "PS05_EVAL:" + hashlib.sha256(canonical_bytes(identity)).hexdigest(),
                **{key: row[key] for key in IDENTITY_FIELDS},
                "positive_score": repr(float(row["positive_score"])),
                "negative_score": repr(float(row["negative_score"])),
                "margin": repr(float(row["margin"])),
                "F_credit": repr(float(row["value"])),
            })

def main() -> None:
    ap = argparse.ArgumentParser(description="Score a complete frozen PS05 leaf ledger.")
    ap.add_argument("--contract-root",type=Path,required=True)
    ap.add_argument("--contract-binding",type=Path,default=Path(__file__).with_name("CONTRACT_BINDING.json"))
    ap.add_argument("--planned-population",type=Path,required=True)
    ap.add_argument("--leaf-ledger",type=Path,required=True)
    ap.add_argument("--delta-num",type=float,required=True)
    ap.add_argument("--reuse-ledger",type=Path,required=True)
    ap.add_argument("--individual-scores",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    a = ap.parse_args()
    if a.delta_num != 0.0: raise ValueError("PS05 formal scorer requires frozen --delta-num 0.0")
    binding = validate_contract(a.contract_root, a.contract_binding)
    rows = load_rows(a.leaf_ledger, a.delta_num)
    population = validate_planned_population(a.planned_population, rows)
    binding.update(validate_relation_population(a.contract_root, a.planned_population))
    trace, formal_f = score(rows, population["scoreable_population_identity_sha256"])
    from bootstrap_ps05_hsf_nested import bootstrap, ci95, DRAWS, SEED
    draws = bootstrap(rows, population["scoreable_population_identity_sha256"], Path(__file__))
    ci = ci95(draws)
    with a.reuse_ledger.open(newline="", encoding="utf-8") as handle:
        reuse = {row["endpoint_id"] for row in csv.DictReader(handle, delimiter="\t") if row["high_reuse_gt_100"] == "True"}
    kept = [row for row in rows if row["anchor_id"] not in reuse and row["positive_endpoint_id"] not in reuse and row["negative_endpoint_id"] not in reuse]
    reuse_f = score(kept, population["scoreable_population_identity_sha256"])[1]
    write_individual(a.individual_scores, rows)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    flat = [dict(x,meta=None) for level in LEVELS for x in trace[level]]
    (a.output_dir/'AGGREGATION_TRACE.jsonl').write_text('\n'.join(canonical_json(x) for x in flat)+'\n',encoding='utf-8', newline='\n')
    summary = {"kind":"PS05_FORMAL_RESULT","status":"FORMAL_COMPLETE","formal_result":True,"formal_interface":"F_ONLY","score_transform_version":"PS05_v1.0.2_F_ONLY","measurement_ledger_sha256":digest(a.leaf_ledger),**population,"delta_num":a.delta_num,"leaf_count":len(rows),"coverage":{"planned":len(rows),"scored":len(rows),"fraction":1.0},"level_counts":{level:len(trace[level]) for level in LEVELS},"F":formal_f,"F_95_CI":list(ci),"bootstrap":{"rng":"PCG64","seed":SEED,"draws":DRAWS},"reuse_sensitivity":{"rule":"reuse_count > 100","kept_leaf_count":len(kept),"F":reuse_f},"individual_scores":{"path":str(a.individual_scores),"sha256":digest(a.individual_scores)},"contract_binding":binding}
    (a.output_dir/'SCORE_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8', newline='\n')
if __name__ == '__main__': main()
