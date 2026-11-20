"""Convert external PS02A representations to the frozen leaf-score ledger."""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if not np.isfinite(denominator):
        raise ValueError("non-finite representation norm")
    if denominator == 0:
        return 0.0
    value = float(left @ right / denominator)
    if not np.isfinite(value):
        raise ValueError("non-finite cosine value")
    return value


def load_representations(path: Path, index_path: Path | None) -> dict[str, np.ndarray]:
    if path.suffix == ".npz":
        if index_path is not None:
            raise ValueError("--representation-index is only valid for safetensors input")
        with np.load(path, allow_pickle=False) as archive:
            if set(archive.files) != {"object_ids", "vectors"}:
                raise ValueError("NPZ must contain exactly object_ids and vectors")
            ids = [str(value) for value in archive["object_ids"]]
            raw_matrix = archive["vectors"]
            if raw_matrix.dtype.kind not in {"i", "u", "f"}:
                raise ValueError("NPZ vectors must use a real numeric dtype")
            matrix = np.asarray(raw_matrix, dtype=np.float32)
        if (
            matrix.ndim != 2
            or matrix.shape[1] <= 0
            or matrix.shape[0] != len(ids)
            or len(ids) != len(set(ids))
            or not np.isfinite(matrix).all()
        ):
            raise ValueError("invalid representation matrix or object IDs")
        return {object_id: matrix[row] for row, object_id in enumerate(ids)}

    if path.suffix != ".safetensors" or index_path is None:
        raise ValueError("input must be NPZ or safetensors with --representation-index")
    try:
        from safetensors import safe_open
    except ImportError as exc:
        raise RuntimeError("safetensors input requires the safetensors package") from exc
    with index_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"occurrence_id", "embedding_key"} <= set(rows[0]):
        raise ValueError("representation index must contain occurrence_id and embedding_key")
    ids = [row["occurrence_id"] for row in rows]
    keys = [row["embedding_key"] for row in rows]
    if len(ids) != len(set(ids)) or len(keys) != len(set(keys)):
        raise ValueError("duplicate occurrence_id or embedding_key in representation index")
    vectors: dict[str, np.ndarray] = {}
    dimension: int | None = None
    with safe_open(str(path), framework="np") as archive:
        if set(archive.keys()) != set(keys):
            raise ValueError("safetensors keys do not match representation index")
        for object_id, key in zip(ids, keys, strict=True):
            raw_vector = archive.get_tensor(key)
            if raw_vector.dtype.kind not in {"i", "u", "f"}:
                raise ValueError(f"representation for {object_id} must use a real numeric dtype")
            vector = np.asarray(raw_vector, dtype=np.float32)
            if vector.ndim != 1 or vector.size == 0 or not np.isfinite(vector).all():
                raise ValueError(f"invalid representation for {object_id}")
            if dimension is None:
                dimension = int(vector.size)
            elif vector.size != dimension:
                raise ValueError("inconsistent representation dimensions")
            vectors[object_id] = vector
    return vectors


def read_fasta_ids(path: Path) -> set[str]:
    identifiers: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].split()[0])
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise ValueError("invalid or duplicate FASTA object IDs")
    return set(identifiers)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundles", type=Path, required=True)
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--representation-index", type=Path)
    parser.add_argument("--objects-fasta", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    vectors = load_representations(args.representations, args.representation_index)
    if args.objects_fasta is not None and set(vectors) != read_fasta_ids(args.objects_fasta):
        raise ValueError("representation IDs do not exactly match the released Forward FASTA")
    with gzip.open(args.bundles, "rt", newline="") as handle:
        bundles = list(csv.DictReader(handle, delimiter="\t"))
    needed = {value for row in bundles for value in (row["anchor"], row["positive"], row["negative"])}
    if not needed <= set(vectors):
        raise ValueError("representations are missing scored bundle endpoints")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["bundle_id", "score_anchor_positive", "score_anchor_negative"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in bundles:
            writer.writerow(
                {
                    "bundle_id": row["bundle_id"],
                    "score_anchor_positive": repr(cosine(vectors[row["anchor"]], vectors[row["positive"]])),
                    "score_anchor_negative": repr(cosine(vectors[row["anchor"]], vectors[row["negative"]])),
                }
            )


if __name__ == "__main__":
    main()
