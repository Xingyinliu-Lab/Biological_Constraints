#!/usr/bin/env python3
"""Materialize the frozen RNA01-02 Forward populations as hash-addressed FASTA."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for sequence_sha256, sequence in sorted(records.items()):
            handle.write(f">{sequence_sha256}\n{sequence}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--object-ledger", required=True, type=Path)
    parser.add_argument("--calibration-ladders", required=True, type=Path)
    parser.add_argument("--l16k-fasta", required=True, type=Path)
    parser.add_argument("--calibration-fasta", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()

    objects = json.loads(args.object_ledger.read_text(encoding="utf-8"))
    if not isinstance(objects, list) or len(objects) != 49152:
        raise ValueError("expected exactly 49,152 formal objects")
    sequences: dict[str, str] = {}
    for row in objects:
        sequence_sha256, sequence = row["sequence_sha256"], row["sequence"]
        if not sequence or set(sequence) - set("ACGU") or sha256_text(sequence) != sequence_sha256:
            raise ValueError(f"invalid sequence for {row.get('object_id', '<unknown>')}")
        existing = sequences.setdefault(sequence_sha256, sequence)
        if existing != sequence:
            raise ValueError(f"conflicting sequence for {sequence_sha256}")
    if len(sequences) != 49149:
        raise ValueError("expected exactly 49,149 unique formal sequences")

    ladders = json.loads(args.calibration_ladders.read_text(encoding="utf-8"))
    calibration_ids = {
        row[key]
        for row in ladders
        for key in ("x0_sequence_sha256", "x1_sequence_sha256", "x2_sequence_sha256")
    }
    if len(ladders) != 128 or len(calibration_ids) != 384 or not calibration_ids <= sequences.keys():
        raise ValueError("calibration panel is not the frozen 128-ladder / 384-endpoint population")
    calibration = {key: sequences[key] for key in calibration_ids}
    write_fasta(args.l16k_fasta, sequences)
    write_fasta(args.calibration_fasta, calibration)
    manifest = {
        "constraint_id": "RNA01-02",
        "formal_object_count": len(objects),
        "l16k_unique_sequence_count": len(sequences),
        "calibration_ladder_count": len(ladders),
        "calibration_unique_sequence_count": len(calibration),
        "l16k_fasta_sha256": hashlib.sha256(args.l16k_fasta.read_bytes()).hexdigest(),
        "calibration_fasta_sha256": hashlib.sha256(args.calibration_fasta.read_bytes()).hexdigest(),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
