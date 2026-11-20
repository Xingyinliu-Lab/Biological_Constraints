#!/usr/bin/env python3
"""Single public PS02A entry from model representations to formal outputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
CONSTRAINT = HERE.parent
REFSET = CONSTRAINT / "refset"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--representation-index", type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    leaf = args.output_dir / "RAW_MEASUREMENTS.tsv.gz"
    individual = args.output_dir / "INDIVIDUAL_SCORES.tsv.gz"
    aggregation = args.output_dir / "DEPENDENCY_AGGREGATION_LEDGER.tsv.gz"
    summary = args.output_dir / "SCORE_SUMMARY.json"
    converter = [
        sys.executable,
        str(HERE / "representations_to_leaf.py"),
        "--bundles",
        str(REFSET / "FORMAL_BUNDLE_LEDGER.tsv.gz"),
        "--representations",
        str(args.representations),
        "--objects-fasta",
        str(REFSET / "endpoint_sequences.faa"),
        "--output",
        str(leaf),
    ]
    if args.representation_index is not None:
        converter.extend(["--representation-index", str(args.representation_index)])
    subprocess.run(converter, check=True)
    subprocess.run(
        [
            sys.executable,
            str(HERE / "score_ps02a_formal.py"),
            "--bundles",
            str(REFSET / "FORMAL_BUNDLE_LEDGER.tsv.gz"),
            "--point-membership",
            str(REFSET / "POINT_DEPENDENCY_MEMBERSHIP.tsv.gz"),
            "--leaf-scores",
            str(leaf),
            "--uncertainty-membership",
            str(REFSET / "UNCERTAINTY_ANCHOR_CLAN_MEMBERSHIP.tsv.gz"),
            "--output",
            str(summary),
            "--individual-scores",
            str(individual),
            "--aggregation-ledger",
            str(aggregation),
        ],
        check=True,
    )
    payload = json.loads(summary.read_text(encoding="utf-8"))
    payload["model_id"] = args.model_id
    payload["formal_interface"] = "F_ONLY"
    payload["score_transform_version"] = "PS02A_v1.0.2_F_ONLY"
    payload["representation_sha256"] = sha256(args.representations)
    payload["representation_index_sha256"] = (
        sha256(args.representation_index) if args.representation_index is not None else None
    )
    payload["refset_assets"] = {
        name: sha256(REFSET / name)
        for name in (
            "FORMAL_BUNDLE_LEDGER.tsv.gz",
            "POINT_DEPENDENCY_MEMBERSHIP.tsv.gz",
            "UNCERTAINTY_ANCHOR_CLAN_MEMBERSHIP.tsv.gz",
        )
    }
    summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    manifest = {"status": "PASS", "files": [
        {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in (leaf, individual, aggregation, summary)
    ]}
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
