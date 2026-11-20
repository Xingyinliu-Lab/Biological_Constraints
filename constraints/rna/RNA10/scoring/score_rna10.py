"""Run the public frozen RNA10 scoring protocol on a compatible measurement ledger."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

from rna10_scoring_core import evaluate, load


def file_record(path: Path) -> dict[str, object]:
    return {
        "path": path.name,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parent
    parser.add_argument("--binding", type=Path, default=root / "RNA10_FORMAL_SCORING_BINDING.json")
    parser.add_argument("--run-plan", type=Path, default=root / "RNA10_FORMAL_RUNPLAN.json")
    parser.add_argument("--measurement-ledger", required=True, type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    plan = json.loads(args.run_plan.read_text(encoding="utf-8"))
    supplied = load(args.measurement_ledger)
    if supplied.get("model_id") not in (None, args.model_id):
        raise ValueError("measurement model identity conflicts with --model-id")
    if not isinstance(supplied.get("records"), list):
        raise ValueError("measurement ledger records missing")
    if not isinstance(supplied.get("representation_dimension"), int) or supplied["representation_dimension"] <= 0:
        raise ValueError("invalid representation dimension")
    if supplied.get("delta_definition") != "z_mut-z_WT":
        raise ValueError("measurement delta definition mismatch")
    if supplied.get("label_encoding") != {"PRESERVING": 0, "DISRUPTING": 1}:
        raise ValueError("measurement label encoding mismatch")

    supplied_lineage = {
        "runplan_sha256": supplied.get("runplan_sha256"),
        "population_identity_sha256": supplied.get("population_identity_sha256"),
        "numeric_method_sha256": supplied.get("numeric_method_sha256"),
        "split_release_sha256": supplied.get("split_release_sha256"),
    }
    source_runplan_sha256 = hashlib.sha256(
        json.dumps(plan, sort_keys=True).encode()
    ).hexdigest()
    expected_supplied_lineage = {
        "runplan_sha256": source_runplan_sha256,
        "population_identity_sha256": binding["population_identity_sha256"],
        "numeric_method_sha256": binding["numeric_method_sha256"],
        "split_release_sha256": binding["split_release_sha256"],
    }
    for key, expected in expected_supplied_lineage.items():
        if supplied_lineage[key] != expected:
            raise ValueError(f"measurement source lineage conflicts with the release: {key}")

    evaluation_plan = copy.deepcopy(plan)
    evaluation_plan["model_interface_id"] = args.model_id
    effective_runplan_sha256 = hashlib.sha256(
        json.dumps(evaluation_plan, sort_keys=True).encode()
    ).hexdigest()
    ledger = copy.deepcopy(supplied)
    ledger["model_id"] = args.model_id
    ledger["runplan_sha256"] = effective_runplan_sha256
    result = evaluate(evaluation_plan, ledger, binding, include_trace=True)
    trace = result.pop("_aggregation_trace")
    result.update(
        {
            "status": "PASS",
            "evaluation_class": "FORMAL_PROTOCOL_EVALUATION",
            "formal_scoring_protocol": True,
            "registered_paper_result": False,
            "constraint_id": "C-RNA-10",
            "model_id": args.model_id,
            "measurement_ledger_sha256": hashlib.sha256(
                args.measurement_ledger.read_bytes()
            ).hexdigest(),
            "supplied_measurement_lineage": supplied_lineage,
            "effective_evaluation_lineage": {
                "runplan_sha256": effective_runplan_sha256,
                "population_identity_sha256": binding["population_identity_sha256"],
                "numeric_method_sha256": binding["numeric_method_sha256"],
                "split_release_sha256": binding["split_release_sha256"],
            },
        }
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    raw_by_leaf = {record["leaf_id"]: record for record in supplied["records"]}
    individual_rows = []
    for scored in trace["leaf_credits"]:
        raw = raw_by_leaf[scored["leaf_id"]]
        individual_rows.append(
            {
                "leaf_id": scored["leaf_id"],
                "matched_block_id": scored["matched_block_id"],
                "source_dependency_unit_id": scored["source_dependency_unit_id"],
                "uncertainty_dependency_cluster_id": scored[
                    "uncertainty_dependency_cluster_id"
                ],
                "h_positive": raw["h_positive"],
                "h_negative": raw["h_negative"],
                "F_credit": scored["F_credit"],
            }
        )
    individual_path = args.out.parent / "INDIVIDUAL_SCORES.json"
    individual_path.write_text(
        json.dumps(individual_rows, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    aggregation_path = args.out.parent / "AGGREGATION_SCORES.json"
    aggregation_path.write_text(
        json.dumps(
            {
                "block_F": trace["block_F"],
                "source_F": trace["source_F"],
                "point_F": trace["point_F"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    bootstrap_path = args.out.parent / "F_BOOTSTRAP_DRAWS.json"
    bootstrap_path.write_text(json.dumps(trace["F_bootstrap"]) + "\n", encoding="utf-8")
    manifest_path = args.out.parent / "OUTPUT_MANIFEST.json"
    manifest_path.write_text(
        json.dumps(
            {
                "constraint_id": "C-RNA-10",
                "model_id": args.model_id,
                "public_entry": Path(__file__).name,
                "outputs": [
                    file_record(args.out),
                    file_record(individual_path),
                    file_record(aggregation_path),
                    file_record(bootstrap_path),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
