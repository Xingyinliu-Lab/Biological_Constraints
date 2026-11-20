#!/usr/bin/env python3
"""Export a PS08 paper-case safetensors/index pair to the public NPZ contract."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from safetensors import safe_open


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_object_ids(population_path: Path) -> set[str]:
    population = json.loads(population_path.read_text(encoding="utf-8-sig"))
    members = population.get("canonical_members")
    if not isinstance(members, list) or len(members) != 4096:
        raise ValueError("paper-case export requires the 4,096-leaf PS08 release population")
    return {
        row[field]
        for row in members
        for field in ("x0_object_id", "x1_object_id", "x2_object_id")
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Export PS08 safetensors representations to object_ids + embeddings NPZ.")
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--population", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.index.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"endpoint_id", "embedding_key"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("embedding index must contain endpoint_id and embedding_key")
        rows = list(reader)
    object_ids = [row["endpoint_id"] for row in rows]
    keys = [row["embedding_key"] for row in rows]
    if not rows or len(object_ids) != len(set(object_ids)) or len(keys) != len(set(keys)):
        raise ValueError("embedding index is empty or contains duplicate IDs/keys")
    expected = expected_object_ids(args.population)
    if set(object_ids) != expected:
        raise ValueError(f"paper-case object universe mismatch: missing={len(expected-set(object_ids))}, extra={len(set(object_ids)-expected)}")

    order = sorted(range(len(object_ids)), key=object_ids.__getitem__)
    ordered_ids = [object_ids[index] for index in order]
    ordered_keys = [keys[index] for index in order]
    vectors: list[np.ndarray] = []
    with safe_open(str(args.embeddings), framework="np") as store:
        if set(store.keys()) != set(keys):
            raise ValueError("safetensors key set does not exactly match the embedding index")
        for key in ordered_keys:
            vector = np.asarray(store.get_tensor(key), dtype=np.float32)
            if vector.ndim != 1 or vector.size == 0 or not np.isfinite(vector).all():
                raise ValueError(f"invalid embedding vector for key {key}")
            vectors.append(vector)
    dimensions = {vector.shape[0] for vector in vectors}
    if len(dimensions) != 1:
        raise ValueError("embedding dimensions are inconsistent")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    matrix = np.stack(vectors)
    np.savez_compressed(args.output, object_ids=np.asarray(ordered_ids, dtype="U128"), embeddings=matrix)
    report = {
        "kind": "PS08_PAPER_CASE_FORWARD_OUTPUT_EXPORT",
        "status": "PASS",
        "row_count": len(ordered_ids),
        "embedding_dimension": matrix.shape[1],
        "safetensors_sha256": sha256(args.embeddings),
        "embedding_index_sha256": sha256(args.index),
        "population_sha256": sha256(args.population),
        "output_sha256": sha256(args.output),
        "output_contract": {"object_ids": list(np.asarray(ordered_ids).shape), "embeddings": list(matrix.shape)},
    }
    args.output.with_name(args.output.name + ".export.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
