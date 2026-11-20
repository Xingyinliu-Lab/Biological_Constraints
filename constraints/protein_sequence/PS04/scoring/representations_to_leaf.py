#!/usr/bin/env python3
"""Convert PS04 whole-domain representations to the frozen leaf ledger."""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np


def as_real_matrix(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values)
    if not np.issubdtype(array.dtype, np.number) or np.iscomplexobj(array):
        raise ValueError("representations must be real numeric vectors")
    return np.asarray(array, dtype=np.float64)


def read_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("released FASTA must contain unique object IDs")
    return ids


def read_index(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"domain_id", "embedding_key"}.issubset(reader.fieldnames):
            raise ValueError("representation index requires domain_id and embedding_key")
        index: dict[str, str] = {}
        keys: set[str] = set()
        for row in reader:
            domain_id, key = row["domain_id"], row["embedding_key"]
            if domain_id in index or key in keys:
                raise ValueError("duplicate domain_id or embedding_key in representation index")
            index[domain_id] = key
            keys.add(key)
    return index


def load_representations(path: Path, index_path: Path | None) -> tuple[list[str], np.ndarray]:
    if path.suffix == ".npz":
        if index_path is not None:
            raise ValueError("--representation-index is only used with safetensors")
        with np.load(path, allow_pickle=False) as archive:
            if set(archive.files) != {"object_ids", "vectors"}:
                raise ValueError("NPZ must contain exactly object_ids and vectors")
            ids = [str(value) for value in archive["object_ids"]]
            matrix = as_real_matrix(archive["vectors"])
    else:
        if index_path is None:
            raise ValueError("safetensors input requires --representation-index")
        try:
            from safetensors import safe_open
        except ImportError as exc:
            raise RuntimeError("install safetensors to read this representation container") from exc
        index = read_index(index_path)
        ids = list(index)
        with safe_open(str(path), framework="np") as store:
            stored_keys = set(store.keys())
            if stored_keys != set(index.values()):
                raise ValueError("safetensors keys and representation index differ")
            matrix = as_real_matrix(np.stack([store.get_tensor(index[value]) for value in ids]))
    if len(ids) != len(set(ids)) or matrix.ndim != 2 or matrix.shape != (len(ids), matrix.shape[1]) or matrix.shape[1] <= 0:
        raise ValueError("representations require unique IDs and a non-empty two-dimensional matrix")
    if not np.isfinite(matrix).all():
        raise ValueError("representations contain NaN or Inf")
    return ids, matrix


def normalize_total(matrix: np.ndarray) -> np.ndarray:
    """Executed PS04 whole-domain cosine: float64 norm with total zero handling."""
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.maximum(norms, 1e-300)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", type=Path, required=True)
    parser.add_argument("--objects-fasta", type=Path, required=True)
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--representation-index", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    ids, matrix = load_representations(args.representations, args.representation_index)
    expected = read_fasta_ids(args.objects_fasta)
    if set(ids) != set(expected) or len(ids) != len(expected):
        raise ValueError("representation object IDs must exactly cover the released Forward FASTA")
    matrix = normalize_total(matrix)
    vectors = {object_id: matrix[index] for index, object_id in enumerate(ids)}

    fields = [
        "negative_endpoint_id", "negative_hsf_id", "negative_topology_id",
        "positive_endpoint_id", "positive_hsf_id", "anchor_id", "anchor_parent_id",
        "anchor_source_unit_id", "anchor_hsf_id", "anchor_topology_id",
        "anchor_architecture_id", "positive_score", "negative_score",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_open = gzip.open if args.output.suffix == ".gz" else open
    count = 0
    with gzip.open(args.relations, "rt", newline="", encoding="utf-8") as source, output_open(
        args.output, "wt", newline="", encoding="utf-8"
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for relation in reader:
            needed = (relation["anchor_id"], relation["positive_id"], relation["negative_id"])
            if any(value not in vectors for value in needed):
                raise ValueError(f"relation references an object absent from the released Forward set: {needed}")
            anchor = vectors[relation["anchor_id"]]
            writer.writerow({
                "negative_endpoint_id": relation["negative_id"],
                "negative_hsf_id": relation["negative_hsf_id"],
                "negative_topology_id": relation["negative_topology_id"],
                "positive_endpoint_id": relation["positive_id"],
                "positive_hsf_id": relation["positive_hsf_id"],
                "anchor_id": relation["anchor_id"],
                "anchor_parent_id": relation["anchor_parent_proxy"],
                "anchor_source_unit_id": relation["anchor_source_unit_id"],
                "anchor_hsf_id": relation["anchor_hsf_id"],
                "anchor_topology_id": relation["anchor_topology_id"],
                "anchor_architecture_id": relation["anchor_architecture_id"],
                "positive_score": repr(float(anchor @ vectors[relation["positive_id"]])),
                "negative_score": repr(float(anchor @ vectors[relation["negative_id"]])),
            })
            count += 1
    if count == 0:
        raise ValueError("formal relation ledger is empty")


if __name__ == "__main__":
    main()
