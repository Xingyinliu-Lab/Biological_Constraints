#!/usr/bin/env python3
"""Single public entry for the frozen PS10 primary evaluation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from ps10_scoring_core import evaluate


SCORING_DIR = Path(__file__).resolve().parent
CONSTRAINT_DIR = SCORING_DIR.parent
REFSET_DIR = CONSTRAINT_DIR / "refset"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def file_record(path: Path) -> dict[str, object]:
    return {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the frozen PS10 primary 22-fold evaluation."
    )
    parser.add_argument("--embeddings", required=True, type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--refset", type=Path, default=REFSET_DIR)
    args = parser.parse_args()

    records = args.refset / "PEPTIDE_REFSET_RECORDS.tsv"
    assignments = args.refset / "PRIMARY_ASSIGNMENTS.tsv"
    method_path = SCORING_DIR / "PS10_NUMERIC_METHOD.json"
    for path in (records, assignments, method_path, args.embeddings):
        if not path.is_file():
            raise FileNotFoundError(path)
    method = json.loads(method_path.read_text(encoding="utf-8"))
    if method.get("status") != "FROZEN_NUMERIC_METHOD":
        raise ValueError("unexpected Numeric Method status")

    calculation = evaluate(records, assignments, args.embeddings)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    individual_path = args.output_dir / "INDIVIDUAL_SCORES.tsv"
    fold_path = args.output_dir / "FOLD_SCORES.tsv"
    summary_path = args.output_dir / "SCORE_SUMMARY.json"
    bootstrap_path = args.output_dir / "BOOTSTRAP_DRAWS.float64le.bin"
    write_tsv(
        individual_path,
        [
            "fold_id",
            "record_id",
            "experimental_classification",
            "binary_label",
            "source_study_units",
            "heldout_decision_score",
        ],
        calculation["individual_rows"],
    )
    write_tsv(
        fold_path,
        ["fold_id", "train_records", "test_records", "purged_records", "heldout_auc", "fold_F"],
        calculation["fold_rows"],
    )
    calculation["bootstrap_draws"].astype("<f8").tofile(bootstrap_path)
    result = {
        "status": "PASS",
        "evaluation_class": "FORMAL_PRIMARY_PROTOCOL_EVALUATION",
        "formal_scoring_protocol": True,
        "registered_paper_result": False,
        "constraint_id": "PS10",
        "model_id": args.model_id,
        "formal_interface": "F_ONLY",
        "score_transform_version": "PS10_v0.1.2_F_ONLY",
        **calculation["summary"],
        "bootstrap": {
            "unit": "heldout_study_fold_F",
            "draws": 20_000,
            "bit_generator": "PCG64",
            "seed": 20_261_008,
            "quantile_method": "linear",
        },
        "input_identity": {
            "records_sha256": sha256(records),
            "primary_assignments_sha256": sha256(assignments),
            "numeric_method_sha256": sha256(method_path),
            "embeddings_sha256": sha256(args.embeddings),
        },
    }
    summary_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    outputs = [individual_path, fold_path, summary_path, bootstrap_path]
    manifest = {
        "constraint_id": "PS10",
        "model_id": args.model_id,
        "formal_interface": "F_ONLY",
        "score_transform_version": "PS10_v0.1.2_F_ONLY",
        "public_entry": Path(__file__).name,
        "outputs": [file_record(path) for path in outputs],
    }
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
