#!/usr/bin/env python3
"""Single public PS04 entry from model representations to formal outputs."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
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

    raw = args.output_dir / "RAW_MEASUREMENTS.tsv"
    individual = args.output_dir / "INDIVIDUAL_SCORES.tsv.gz"
    convert = [
        sys.executable, str(HERE / "representations_to_leaf.py"),
        "--relations", str(REFSET / "FORMAL_RELATION_LEDGER.tsv.gz"),
        "--objects-fasta", str(REFSET / "FORMAL_ENDPOINTS.faa"),
        "--representations", str(args.representations),
        "--output", str(raw),
    ]
    if args.representation_index is not None:
        convert.extend(["--representation-index", str(args.representation_index)])
    subprocess.run(convert, check=True)
    subprocess.run([
        sys.executable, str(HERE / "score_ps04_formal.py"),
        "--binding", str(HERE / "FORMAL_BINDING.json"),
        "--formal-relations", str(REFSET / "FORMAL_RELATION_LEDGER.tsv.gz"),
        "--leaf-ledger", str(raw),
        "--uncertainty-membership", str(REFSET / "UNCERTAINTY_ANCHOR_SOURCE_UNIT_MEMBERSHIP.tsv.gz"),
        "--endpoint-reuse-ledger", str(REFSET / "ENDPOINT_REUSE_LEDGER.tsv.gz"),
        "--individual-scores", str(individual),
        "--output-dir", str(args.output_dir),
    ], check=True)

    summary = args.output_dir / "FORMAL_RESULT.json"
    payload = json.loads(summary.read_text(encoding="utf-8"))
    payload["model_id"] = args.model_id
    payload["formal_interface"] = "F_ONLY"
    payload["score_transform_version"] = "PS04_v1.0.2_F_ONLY"
    payload["representation_sha256"] = sha256(args.representations)
    payload["representation_index_sha256"] = sha256(args.representation_index) if args.representation_index else None
    payload["formal_result"] = True
    summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    raw_compressed = args.output_dir / "RAW_MEASUREMENTS.tsv.gz"
    with raw.open("rb") as source, gzip.open(raw_compressed, "wb", compresslevel=6) as destination:
        shutil.copyfileobj(source, destination, length=1 << 20)
    raw.unlink()
    aggregation_trace = args.output_dir / "AGGREGATION_TRACE.jsonl.gz"
    outputs = (raw_compressed, individual, aggregation_trace, summary)
    manifest = {"status": "PASS", "files": [
        {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs
    ]}
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
