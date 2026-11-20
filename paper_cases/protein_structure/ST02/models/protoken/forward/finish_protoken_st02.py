#!/usr/bin/env python3
"""Convert the registered ProToken codes to ST02 scorer input."""

from __future__ import annotations

import argparse
import csv
import json
import pickle
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tokens", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--objects", type=Path, required=True)
    parser.add_argument("--relations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.tokens.open("rb") as handle:
        raw_tokens = pickle.load(handle)
    with args.objects.open(encoding="utf-8-sig", newline="") as handle:
        objects = list(csv.DictReader(handle, delimiter="\t"))
    by_locator_stem = {
        Path(row["coordinate_locator"]).stem: row["record_id"] for row in objects
    }
    formal_ids = set(by_locator_stem.values())
    tokens: dict[str, np.ndarray] = {}
    for key, value in raw_tokens.items():
        raw_id = str(key).removesuffix("_clean")
        object_id = raw_id if raw_id in formal_ids else by_locator_stem.get(raw_id)
        if object_id is None or object_id in tokens:
            raise ValueError("ProToken token/object mapping failure")
        tokens[object_id] = np.asarray(value, dtype=np.int64)
    with args.checkpoint.open("rb") as handle:
        checkpoint = pickle.load(handle)
    codebook = np.asarray(
        checkpoint["params"]["vq_tokenizer"]["codebook"], dtype=np.float32
    )
    with args.relations.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        object_ids = sorted({
            row[key]
            for row in rows
            for key in ("anchor_object_id", "true_target_object_id", "negative_object_id")
        })
    if set(object_ids) != set(tokens):
        raise ValueError("ProToken object closure failure")
    vectors = np.asarray(
        [codebook[tokens[object_id]].mean(axis=0) for object_id in object_ids],
        dtype=np.float32,
    )
    if vectors.ndim != 2 or vectors.shape[1] == 0 or not np.isfinite(vectors).all():
        raise ValueError("invalid ProToken embedding matrix")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        object_ids=np.asarray(object_ids),
        embeddings=vectors,
    )
    print(json.dumps({
        "status": "PASS",
        "objects": len(object_ids),
        "dimension": int(vectors.shape[1]),
        "pooling": "unweighted mean of registered ProToken codebook vectors",
        "output": str(args.output),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
