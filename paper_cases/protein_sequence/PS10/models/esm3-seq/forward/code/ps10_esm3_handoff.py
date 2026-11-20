#!/usr/bin/env python3
"""Prepare PS10 ESM3 FASTA input and export its tensors to the public NPZ contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_records(path: Path) -> list[dict[str, str]]:
    git_root = Path(__file__).resolve().parents[7]
    scoring_dir = git_root / "constraints" / "protein_sequence" / "PS10" / "scoring"
    sys.path.insert(0, str(scoring_dir))
    from ps10_refset import validate_records

    return validate_records(path)


def prepare_fasta(records_path: Path, output: Path) -> None:
    records = load_records(records_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="ascii", newline="\n") as handle:
        for row in records:
            handle.write(f">{row['record_id']}\n{row['normalized_peptide_identity']}\n")
    report = {
        "kind": "PS10_ESM3_FASTA_EXPORT",
        "status": "PASS",
        "record_count": len(records),
        "records_sha256": sha256(records_path),
        "fasta_sha256": sha256(output),
        "sequence_field": "normalized_peptide_identity",
    }
    output.with_name(output.name + ".export.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )


def export_npz(records_path: Path, embeddings_path: Path, index_path: Path, output: Path) -> None:
    try:
        from safetensors import safe_open
    except ImportError as exc:
        raise RuntimeError("export-npz requires the safetensors package") from exc

    records = load_records(records_path)
    record_by_id = {row["record_id"]: row for row in records}
    with index_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "record_id",
            "embedding_key",
            "sequence_length",
            "embedding_dim",
            "embedding_sha256",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("invalid ESM3 embedding index")
        index_rows = list(reader)
    ids = [row["record_id"] for row in index_rows]
    keys = [row["embedding_key"] for row in index_rows]
    if len(ids) != len(set(ids)) or len(keys) != len(set(keys)):
        raise ValueError("duplicate ESM3 record ID or tensor key")
    if set(ids) != set(record_by_id):
        raise ValueError("ESM3 index/PS10 record population mismatch")

    vectors: dict[str, np.ndarray] = {}
    with safe_open(str(embeddings_path), framework="np") as store:
        if set(store.keys()) != set(keys):
            raise ValueError("ESM3 tensor/index key closure failure")
        for row in index_rows:
            identity = row["record_id"]
            if int(row["sequence_length"]) != int(record_by_id[identity]["length"]):
                raise ValueError(f"ESM3 sequence length mismatch: {identity}")
            vector = np.asarray(store.get_tensor(row["embedding_key"]), dtype=np.float32)
            if vector.ndim != 1 or vector.size == 0 or not np.isfinite(vector).all():
                raise ValueError(f"invalid ESM3 vector: {identity}")
            if vector.size != int(row["embedding_dim"]):
                raise ValueError(f"ESM3 embedding dimension mismatch: {identity}")
            if hashlib.sha256(np.ascontiguousarray(vector).tobytes()).hexdigest() != row[
                "embedding_sha256"
            ]:
                raise ValueError(f"ESM3 vector hash mismatch: {identity}")
            vectors[identity] = vector
    dimensions = {vector.size for vector in vectors.values()}
    if len(dimensions) != 1:
        raise ValueError("ESM3 embedding dimensions are inconsistent")
    ordered_ids = [row["record_id"] for row in records]
    matrix = np.stack([vectors[identity] for identity in ordered_ids])
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        record_id=np.asarray(ordered_ids, dtype="U64"),
        embedding=matrix,
    )
    report = {
        "kind": "PS10_ESM3_PUBLIC_NPZ_EXPORT",
        "status": "PASS",
        "record_count": len(ordered_ids),
        "embedding_dimension": matrix.shape[1],
        "records_sha256": sha256(records_path),
        "safetensors_sha256": sha256(embeddings_path),
        "embedding_index_sha256": sha256(index_path),
        "output_sha256": sha256(output),
    }
    output.with_name(output.name + ".export.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    fasta_parser = subparsers.add_parser("prepare-fasta")
    fasta_parser.add_argument("--records", required=True, type=Path)
    fasta_parser.add_argument("--output", required=True, type=Path)
    export_parser = subparsers.add_parser("export-npz")
    export_parser.add_argument("--records", required=True, type=Path)
    export_parser.add_argument("--embeddings", required=True, type=Path)
    export_parser.add_argument("--index", required=True, type=Path)
    export_parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.action == "prepare-fasta":
        prepare_fasta(args.records, args.output)
    else:
        export_npz(args.records, args.embeddings, args.index, args.output)


if __name__ == "__main__":
    main()
