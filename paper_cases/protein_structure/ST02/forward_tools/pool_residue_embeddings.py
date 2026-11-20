#!/usr/bin/env python3
"""Pool registered per-residue structure representations to ST02 object vectors."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def formal_object_ids(relations: Path) -> list[str]:
    with relations.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return sorted({
            row[key]
            for row in rows
            for key in ("anchor_object_id", "true_target_object_id", "negative_object_id")
        })


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply the registered unweighted residue mean and emit ST02 scorer input."
    )
    parser.add_argument("--relations", type=Path, required=True)
    parser.add_argument("--objects", type=Path, required=True)
    parser.add_argument("--residue-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    object_ids = formal_object_ids(args.relations)
    with args.objects.open(encoding="utf-8-sig", newline="") as handle:
        locator = {row["record_id"]: Path(row["coordinate_locator"]).stem for row in csv.DictReader(handle, delimiter="\t")}
    vectors = []
    dimension = None
    for object_id in object_ids:
        if object_id not in locator:
            raise ValueError(f"formal object missing from objects table: {object_id}")
        path = args.residue_dir / f"{locator[object_id]}.npz"
        if not path.is_file():
            raise FileNotFoundError(path)
        with np.load(path, allow_pickle=False) as saved:
            residues = np.asarray(saved["residue_embedding"])
        if residues.ndim != 2 or residues.shape[0] == 0 or residues.shape[1] == 0:
            raise ValueError(f"invalid residue embedding shape: {path}")
        if not np.isfinite(residues).all():
            raise ValueError(f"non-finite residue embedding: {path}")
        if dimension is None:
            dimension = int(residues.shape[1])
        elif residues.shape[1] != dimension:
            raise ValueError(f"embedding dimension mismatch: {path}")
        vectors.append(residues.mean(axis=0, dtype=np.float64).astype(np.float32))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        object_ids=np.asarray(object_ids),
        embeddings=np.asarray(vectors, dtype=np.float32),
    )
    print(json.dumps({
        "status": "PASS",
        "objects": len(object_ids),
        "dimension": dimension,
        "pooling": "unweighted residue mean in float64, output float32",
        "output": str(args.output),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
