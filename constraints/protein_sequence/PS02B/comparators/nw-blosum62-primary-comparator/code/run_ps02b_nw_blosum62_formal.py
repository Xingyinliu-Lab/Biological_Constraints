#!/usr/bin/env python3
"""Run frozen PS02B NW-BLOSUM62 measurements through the public score core."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
COMPARATOR = HERE.parent
CONSTRAINT = COMPARATOR.parents[1]
REFSET = CONSTRAINT / "refset"
CONTRACT = COMPARATOR / "method" / "PS02B_NW_BLOSUM62_COMPARATOR_CONTRACT_v0.1.0_FROZEN.json"
MEASURER = HERE / "ps02b_nw_blosum62_comparator.py"
SCORER = CONSTRAINT / "scoring" / "score_ps02b_formal.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256(path)
    if actual != expected:
        raise SystemExit(f"frozen binding mismatch for {label}: {actual}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--needle", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"output directory already exists: {args.output_dir}")

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    frozen = contract["frozen_execution_bindings"]
    relations = REFSET / "FORMAL_RELATION_LEDGER.tsv.gz"
    objects = REFSET / "FORMAL_OBJECT_LEDGER.tsv.gz"
    fasta = REFSET / "endpoint_sequences.faa"
    points = REFSET / "POINT_DEPENDENCY_MEMBERSHIP.tsv.gz"
    uncertainty = REFSET / "UNCERTAINTY_GROUP_MEMBERSHIP.tsv.gz"
    for path, key, label in (
        (relations, "formal_relation_ledger_sha256", "relation ledger"),
        (objects, "formal_object_ledger_sha256", "object ledger"),
        (fasta, "endpoint_fasta_sha256", "endpoint FASTA"),
        (points, "point_membership_sha256", "point membership"),
        (uncertainty, "uncertainty_membership_sha256", "uncertainty membership"),
        (MEASURER, "comparator_code_sha256", "comparator implementation"),
    ):
        require_hash(path, frozen[key], label)

    args.output_dir.mkdir(parents=True)
    measurement = args.output_dir / "measurement"
    subprocess.run(
        [sys.executable, str(MEASURER), "--bundles", str(relations), "--fasta", str(fasta), "--needle", str(args.needle), "--outdir", str(measurement), "--workers", str(args.workers)],
        check=True,
    )
    result = args.output_dir / "SCORE_SUMMARY.json"
    individual = args.output_dir / "INDIVIDUAL_SCORES.tsv.gz"
    subprocess.run(
        [sys.executable, str(SCORER), "--relations", str(relations), "--objects", str(objects), "--point-membership", str(points), "--uncertainty-membership", str(uncertainty), "--leaf-scores", str(measurement / "FORMAL_LEAF_SCORES.tsv.gz"), "--output", str(result), "--individual-scores", str(individual)],
        check=True,
    )
    release = {
        "artifact_type": "AuxiliaryMethodExecutionRelease",
        "status": "FORMAL_EXECUTION_COMPLETE",
        "comparator_method_id": contract["comparator_method_id"],
        "contract_sha256": sha256(CONTRACT),
        "public_score_core_sha256": sha256(SCORER),
        "emboss_needle_version_required": contract["method"]["tool_version"],
        "historical_needle_sha256": frozen["needle_sha256"],
        "executed_needle_sha256": sha256(args.needle),
        "score_ledger_sha256": sha256(measurement / "FORMAL_LEAF_SCORES.tsv.gz"),
        "individual_scores_sha256": sha256(individual),
        "result_sha256": sha256(result),
        "workers": args.workers,
    }
    (args.output_dir / "EXECUTION_RELEASE.json").write_text(json.dumps(release, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
