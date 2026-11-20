#!/usr/bin/env python3
"""Pool per-coordinate residue arrays into the canonical ST01 NPZ contract."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--objects", required=True, type=Path)
    parser.add_argument("--residue-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.objects.open(encoding="utf-8", newline="") as handle:
        objects = list(csv.DictReader(handle, delimiter="\t"))
    if not objects or set(objects[0]) < {"object_id", "coordinate_locator"}:
        raise ValueError("objects.tsv lacks canonical object mapping")

    ids, vectors = [], []
    for row in objects:
        object_id = row["object_id"]
        stem = Path(row["coordinate_locator"]).stem
        path = args.residue_dir / f"{stem}.npz"
        with np.load(path, allow_pickle=False) as data:
            values = np.asarray(data["residue_embedding"])
        if values.ndim != 2 or values.shape[0] == 0 or not np.isfinite(values).all():
            raise ValueError(f"invalid residue embedding: {path}")
        ids.append(object_id)
        vectors.append(values.mean(axis=0, dtype=np.float64).astype(np.float32))
    matrix = np.asarray(vectors, dtype=np.float32)
    if np.any(np.linalg.norm(matrix, axis=1) == 0):
        raise ValueError("zero-norm pooled representation")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output, object_ids=np.asarray(ids), embeddings=matrix)


if __name__ == "__main__":
    main()
