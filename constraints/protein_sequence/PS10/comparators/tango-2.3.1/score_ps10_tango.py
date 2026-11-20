#!/usr/bin/env python3
"""Apply the shared PS10 fold aggregation to frozen TANGO scores."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score


PS10_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PS10_DIR / "scoring"))
from ps10_scoring_core import summarize_fold_f  # noqa: E402
from ps10_refset import validate_refset  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], values: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(values)


def file_record(path: Path) -> dict[str, object]:
    return {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", required=True, type=Path)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--numeric-method", required=True, type=Path)
    parser.add_argument("--specification", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    closure = validate_refset(args.records, args.ledger)
    records = closure["records"]
    ledger = closure["assignments"]
    score_rows = rows(args.scores)
    specification = json.loads(args.specification.read_text())
    method = json.loads(args.numeric_method.read_text())
    if specification["status"] != "FROZEN_COMPARATOR_SPECIFICATION" or method["status"] != "FROZEN_NUMERIC_METHOD":
        raise RuntimeError("unexpected comparator specification or Numeric Method status")

    actual = {
        "records_sha256": sha256(args.records),
        "ledger_sha256": sha256(args.ledger),
        "scores_sha256": sha256(args.scores),
        "numeric_method_sha256": sha256(args.numeric_method),
    }
    identifiers = {row["record_id"] for row in records}
    scores = {row["record_id"]: float(row["tango_aggregation"]) for row in score_rows}
    if len(scores) != len(identifiers) or set(scores) != identifiers or not np.isfinite(list(scores.values())).all():
        raise RuntimeError("white-box score population mismatch")
    if len(ledger) != len(identifiers) * 22 or len({(row["fold_id"], row["record_id"]) for row in ledger}) != len(ledger):
        raise RuntimeError("ledger invariant failed")

    labels = {
        row["record_id"]: int(row["experimental_classification"] == "AMYLOID")
        for row in records
    }
    fold_rows: list[dict[str, object]] = []
    individual_rows: list[dict[str, object]] = []
    fold_f: list[float] = []
    for fold_id in sorted({row["fold_id"] for row in ledger}):
        test = [
            row["record_id"]
            for row in ledger
            if row["fold_id"] == fold_id and row["role"] == "test"
        ]
        outcomes = np.asarray([labels[key] for key in test], dtype=int)
        if not test or len(np.unique(outcomes)) != 2:
            raise RuntimeError("unevaluable held-out dependency unit")
        auc = float(roc_auc_score(outcomes, [scores[key] for key in test]))
        fold_score = 2.0 * auc - 1.0
        fold_f.append(fold_score)
        fold_rows.append({"fold_id": fold_id, "test_records": len(test), "heldout_auc": repr(auc), "fold_F": repr(fold_score)})
        individual_rows.extend(
            {
                "fold_id": fold_id,
                "record_id": record_id,
                "binary_label": labels[record_id],
                "tango_aggregation": repr(scores[record_id]),
            }
            for record_id in test
        )

    summary, draws = summarize_fold_f(fold_f)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.output_dir / "INDIVIDUAL_SCORES.tsv",
        ["fold_id", "record_id", "binary_label", "tango_aggregation"],
        individual_rows,
    )
    write_tsv(
        args.output_dir / "FOLD_SCORES.tsv",
        ["fold_id", "test_records", "heldout_auc", "fold_F"],
        fold_rows,
    )
    draws.astype("<f8").tofile(args.output_dir / "BOOTSTRAP_DRAWS.float64le.bin")
    output = {
        "status": "PASS",
        "constraint_id": "PS10",
        "comparator_id": "tango-2.3.1",
        "formal_interface": "F_ONLY",
        "score_transform_version": "PS10_v0.1.2_F_ONLY",
        **summary,
        "fold_F": fold_rows,
        "coverage": {
            "records": len(records),
            "folds": len(fold_rows),
            "heldout_prediction_rows": len(individual_rows),
        },
        "refset_closure": closure["counts"],
        **actual,
        "specification_sha256": sha256(args.specification),
    }
    summary_path = args.output_dir / "SCORE_SUMMARY.json"
    summary_path.write_text(
        json.dumps(output, indent=2) + "\n", encoding="utf-8"
    )
    outputs = [
        args.output_dir / "INDIVIDUAL_SCORES.tsv",
        args.output_dir / "FOLD_SCORES.tsv",
        summary_path,
        args.output_dir / "BOOTSTRAP_DRAWS.float64le.bin",
    ]
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(
        json.dumps({"constraint_id": "PS10", "comparator_id": "tango-2.3.1", "formal_interface": "F_ONLY", "score_transform_version": "PS10_v0.1.2_F_ONLY", "outputs": [file_record(path) for path in outputs]}, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
