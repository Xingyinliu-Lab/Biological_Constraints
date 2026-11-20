#!/usr/bin/env python3
"""Single public PS02B entry from model representations to formal outputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REFSET = HERE.parent / "refset"


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
    raw = args.output_dir / "RAW_MEASUREMENTS.tsv.gz"
    individual = args.output_dir / "INDIVIDUAL_SCORES.tsv.gz"
    summary = args.output_dir / "SCORE_SUMMARY.json"
    convert = [
        sys.executable,
        str(HERE / "representations_to_leaf.py"),
        "--relations",
        str(REFSET / "FORMAL_RELATION_LEDGER.tsv.gz"),
        "--embeddings",
        str(args.representations),
        "--objects-fasta",
        str(REFSET / "endpoint_sequences.faa"),
        "--output",
        str(raw),
        "--model-id",
        args.model_id,
    ]
    if args.representation_index is not None:
        convert.extend(["--embedding-index", str(args.representation_index)])
    subprocess.run(convert, check=True)
    subprocess.run(
        [
            sys.executable,
            str(HERE / "score_ps02b_formal.py"),
            "--relations",
            str(REFSET / "FORMAL_RELATION_LEDGER.tsv.gz"),
            "--objects",
            str(REFSET / "FORMAL_OBJECT_LEDGER.tsv.gz"),
            "--point-membership",
            str(REFSET / "POINT_DEPENDENCY_MEMBERSHIP.tsv.gz"),
            "--uncertainty-membership",
            str(REFSET / "UNCERTAINTY_GROUP_MEMBERSHIP.tsv.gz"),
            "--leaf-scores",
            str(raw),
            "--output",
            str(summary),
            "--individual-scores",
            str(individual),
        ],
        check=True,
    )
    payload = json.loads(summary.read_text(encoding="utf-8"))
    payload["model_id"] = args.model_id
    payload["formal_interface"] = "F_ONLY"
    payload["score_transform_version"] = "PS02B_v1.0.2_F_ONLY"
    payload["representation_sha256"] = sha256(args.representations)
    payload["representation_index_sha256"] = sha256(args.representation_index) if args.representation_index else None
    summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    outputs = (raw, individual, summary)
    manifest = {"status": "PASS", "files": [
        {"path": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)} for p in outputs
    ]}
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
