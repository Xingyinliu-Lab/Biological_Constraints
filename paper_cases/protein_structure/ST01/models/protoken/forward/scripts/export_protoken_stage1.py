#!/usr/bin/env python3
"""Export canonical ST01 embeddings from ProToken stage-1 VQ indices."""
from __future__ import annotations

import argparse
import csv
import pickle
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--token-index", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--objects", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.token_index.open("rb") as handle:
        raw_tokens = pickle.load(handle)
    with args.checkpoint.open("rb") as handle:
        codebook = np.asarray(pickle.load(handle)["params"]["vq_tokenizer"]["codebook"], dtype=np.float32)
    with args.objects.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row["object_id"] for row in rows]
    coordinate_to_id = {Path(row["coordinate_locator"]).stem: row["object_id"] for row in rows}
    if len(coordinate_to_id) != len(rows):
        raise ValueError("ST01 coordinate filename stems are not unique")
    tokens = {}
    for key, value in raw_tokens.items():
        coordinate_stem = Path(str(key)).stem.removesuffix("_clean")
        if coordinate_stem not in coordinate_to_id:
            raise ValueError(f"ProToken stage-1 key is not a frozen ST01 coordinate: {key}")
        object_id = coordinate_to_id[coordinate_stem]
        if object_id in tokens:
            raise ValueError(f"duplicate ProToken stage-1 object: {object_id}")
        tokens[object_id] = np.asarray(value, dtype=np.int64)
    if set(ids) != set(tokens):
        raise ValueError("ProToken stage-1 object coverage does not match frozen ST01 objects")
    vectors = np.asarray([codebook[tokens[object_id]].mean(axis=0) for object_id in ids], dtype=np.float32)
    if vectors.ndim != 2 or not np.isfinite(vectors).all() or np.any(np.linalg.norm(vectors, axis=1) == 0):
        raise ValueError("invalid ProToken global representations")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output, object_ids=np.asarray(ids), embeddings=vectors)


if __name__ == "__main__":
    main()
