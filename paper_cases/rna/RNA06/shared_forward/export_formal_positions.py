#!/usr/bin/env python3
"""Export sequence matrices to the frozen RNA06 sparse position contract."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_matrices(path):
    path = Path(path)
    if path.suffix.lower() == ".npz":
        archive = np.load(path, allow_pickle=False)
        matrices = {str(key): np.asarray(archive[key]) for key in archive.files}
    else:
        import torch

        value = torch.load(path, map_location="cpu", weights_only=True)
        if not isinstance(value, dict):
            raise ValueError("matrix input must be a sequence-keyed mapping")
        matrices = {
            str(key): item.detach().cpu().numpy() if isinstance(item, torch.Tensor) else np.asarray(item)
            for key, item in value.items()
        }
    if not matrices:
        raise ValueError("empty sequence-matrix input")
    for key, matrix in matrices.items():
        if matrix.dtype.kind not in "iuf" or matrix.ndim != 2 or matrix.shape[1] <= 0:
            raise ValueError(f"invalid real sequence matrix: {key}")
        if not np.isfinite(matrix).all():
            raise ValueError(f"non-finite sequence matrix: {key}")
    return matrices


def main():
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[4] / "constraints" / "rna" / "RNA06" / "refset"
    parser.add_argument("--matrices", required=True)
    parser.add_argument("--objects", default=str(root / "rna_objects.tsv"))
    parser.add_argument("--directed", default=str(root / "RNA06_FORMAL_DIRECTED_ENDPOINT_LEDGER.tsv"))
    parser.add_argument("--candidates", default=str(root / "RNA06_FORMAL_PAIRED_NONPARTNER_CANDIDATE_LEDGER.tsv"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    matrices = load_matrices(args.matrices)
    sequence_by_record = {row["record_id"]: row["sequence_sha256"] for row in read_tsv(args.objects)}
    required = set()
    for row in read_tsv(args.directed):
        sequence = sequence_by_record[row["record_id"]]
        required.add((sequence, int(row["anchor_position0"])))
        required.add((sequence, int(row["true_partner_position0"])))
    for row in read_tsv(args.candidates):
        required.add((sequence_by_record[row["record_id"]], int(row["paired_nonpartner_position0"])))
    if set(matrices) != set(sequence_by_record.values()):
        raise ValueError("sequence matrices must exactly cover the released unique RNA sequences")
    ordered = sorted(required)
    dimensions = {matrix.shape[1] for matrix in matrices.values()}
    if len(dimensions) != 1:
        raise ValueError("sequence matrices have inconsistent dimensions")
    vectors = []
    for sequence, position in ordered:
        matrix = matrices[sequence]
        if position < 0 or position >= matrix.shape[0]:
            raise ValueError(f"position outside sequence matrix: {sequence}:{position}")
        vectors.append(matrix[position])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        sequence_sha256=np.asarray([key[0] for key in ordered]),
        position0=np.asarray([key[1] for key in ordered], dtype=np.int64),
        vectors=np.asarray(vectors, dtype=np.float32),
    )


if __name__ == "__main__":
    main()
