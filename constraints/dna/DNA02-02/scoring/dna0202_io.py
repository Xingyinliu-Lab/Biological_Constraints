from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    key = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if key is not None:
                if key in records:
                    raise ValueError("duplicate FASTA identifier")
                records[key] = "".join(chunks).upper()
            key, chunks = raw[1:].split()[0], []
        else:
            if key is None:
                raise ValueError("sequence before FASTA header")
            chunks.append(raw.strip())
    if key is not None:
        if key in records:
            raise ValueError("duplicate FASTA identifier")
        records[key] = "".join(chunks).upper()
    return records


def validate_sequence_universe(sequence_objects_path: Path, fasta_path: Path) -> set[str]:
    rows = read_tsv(sequence_objects_path)
    fasta = read_fasta(fasta_path)
    ids: set[str] = set()
    for row in rows:
        object_id = row["object_id"]
        sequence = fasta.get(object_id)
        if object_id in ids or sequence is None:
            raise ValueError("duplicate or missing sequence object")
        if len(sequence) != int(row["sequence_length"]) or hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
            raise ValueError("sequence identity mismatch")
        if set(sequence) - set("ACGT"):
            raise ValueError("non-ACGT sequence")
        ids.add(object_id)
    if ids != set(fasta):
        raise ValueError("FASTA and sequence-object universes differ")
    return ids


def load_embeddings(matrix_path: Path, index_path: Path | None) -> dict[str, np.ndarray]:
    if matrix_path.suffix == ".npz":
        archive = np.load(matrix_path, allow_pickle=False)
        if set(archive.files) != {"object_ids", "vectors"}:
            raise ValueError("NPZ must contain exactly object_ids and vectors")
        object_ids = [str(value) for value in archive["object_ids"]]
        matrix = np.asarray(archive["vectors"])
        if matrix.ndim != 2 or matrix.shape[1] <= 0 or matrix.shape[0] != len(object_ids) or len(set(object_ids)) != len(object_ids):
            raise ValueError("invalid NPZ embedding dimensions or duplicate objects")
        if not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(matrix.dtype, np.complexfloating) or np.issubdtype(matrix.dtype, np.bool_):
            raise ValueError("embeddings must use a real numeric dtype")
        if not np.isfinite(matrix).all():
            raise ValueError("non-finite embedding")
        return {object_id: np.asarray(matrix[index], dtype=np.float64) for index, object_id in enumerate(object_ids)}
    if index_path is None:
        raise ValueError("embedding index is required for NPY input")
    rows = read_tsv(index_path)
    required = {"object_id", "row_index", "sequence_sha256"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("invalid embedding index header")
    matrix = np.load(matrix_path, mmap_mode="r")
    if matrix.ndim != 2 or matrix.shape[1] <= 0 or matrix.shape[0] != len(rows):
        raise ValueError("embedding matrix/index dimension mismatch")
    if not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(matrix.dtype, np.complexfloating) or np.issubdtype(matrix.dtype, np.bool_):
        raise ValueError("embeddings must use a real numeric dtype")
    result: dict[str, np.ndarray] = {}
    seen_rows: set[int] = set()
    for row in rows:
        object_id = row["object_id"]
        index = int(row["row_index"])
        if object_id in result or index in seen_rows or index < 0 or index >= matrix.shape[0]:
            raise ValueError("duplicate or invalid embedding index")
        vector = np.asarray(matrix[index], dtype=np.float64)
        if not np.isfinite(vector).all():
            raise ValueError(f"non-finite embedding: {object_id}")
        result[object_id] = vector
        seen_rows.add(index)
    return result


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
