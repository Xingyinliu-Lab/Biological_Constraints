#!/usr/bin/env python3
"""Convert PS05 whole-domain representations into the frozen leaf ledger."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def as_real_matrix(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values)
    if not np.issubdtype(array.dtype, np.number) or np.iscomplexobj(array):
        raise ValueError("representations must be real numeric vectors")
    return np.asarray(array, dtype=np.float64)

FIELDS = [
    "negative_endpoint_id", "negative_hsf_id", "positive_endpoint_id", "positive_hsf_id",
    "anchor_id", "anchor_parent_id", "anchor_source_unit_id", "anchor_hsf_id",
    "anchor_topology_id", "anchor_architecture_id", "assay_id", "planned_bundle_ref",
    "interface_coordinate_ref", "execution_attempt_ref", "measurement_contract_ref",
    "positive_score", "negative_score",
]
IDENTITY = FIELDS[:10]


def fasta_ids(path: Path) -> list[str]:
    ids = [line[1:].strip().split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">")]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("released FASTA must contain unique object IDs")
    return ids


def load_vectors(path: Path, index_path: Path | None) -> tuple[list[str], np.ndarray]:
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
        with index_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or not {"domain_id", "embedding_key"}.issubset(reader.fieldnames):
                raise ValueError("representation index requires domain_id and embedding_key")
            index: dict[str, str] = {}
            keys: set[str] = set()
            for row in reader:
                if row["domain_id"] in index or row["embedding_key"] in keys:
                    raise ValueError("duplicate domain_id or embedding_key")
                index[row["domain_id"]] = row["embedding_key"]
                keys.add(row["embedding_key"])
        ids = list(index)
        with safe_open(str(path), framework="np") as store:
            if set(store.keys()) != set(index.values()):
                raise ValueError("safetensors keys and index differ")
            matrix = as_real_matrix(np.stack([store.get_tensor(index[value]) for value in ids]))
    if len(ids) != len(set(ids)) or matrix.ndim != 2 or matrix.shape[0] != len(ids) or matrix.shape[1] <= 0 or not np.isfinite(matrix).all():
        raise ValueError("invalid representation matrix or object IDs")
    norms = np.linalg.norm(matrix, axis=1)
    if np.any(norms <= 0):
        raise ValueError("PS05 requires finite nonzero whole-domain vectors")
    return ids, matrix / norms[:, None]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", type=Path, required=True)
    parser.add_argument("--planned-population", type=Path, required=True)
    parser.add_argument("--objects-fasta", type=Path, required=True)
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--representation-index", type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    relations = list(csv.DictReader(args.relations.open(encoding="utf-8", newline=""), delimiter="\t"))
    planned = json.loads(args.planned_population.read_text(encoding="utf-8"))
    planned_identity = {tuple(values[:5]): values for values in planned["leaf_identities"]}
    if len(planned_identity) != len(planned["leaf_identities"]):
        raise ValueError("duplicate identity in planned population")
    object_ids, matrix = load_vectors(args.representations, args.representation_index)
    expected = fasta_ids(args.objects_fasta)
    if set(object_ids) != set(expected) or len(object_ids) != len(expected):
        raise ValueError("representations must exactly cover the released Forward FASTA")
    vectors = {object_id: matrix[index] for index, object_id in enumerate(object_ids)}

    rows = []
    for relation in relations:
        key = (relation["negative_domain_id"], relation["negative_hsf_id"], relation["positive_domain_id"], relation["positive_hsf_id"], relation["anchor_domain_id"])
        if key not in planned_identity:
            raise ValueError("relation absent from planned population")
        identity = planned_identity[key]
        row = dict(zip(IDENTITY, identity))
        anchor = vectors[relation["anchor_domain_id"]]
        row.update({
            "assay_id": "C-PS-05-A",
            "planned_bundle_ref": "RS-PS05-L_4K",
            "interface_coordinate_ref": f"PS05_{args.model_id}_GLOBAL_DOMAIN_SAFE_COSINE_v1",
            "execution_attempt_ref": f"PS05_{args.model_id}_public_evaluation",
            "measurement_contract_ref": "PS05_MEASUREMENT_v0.1.0",
            "positive_score": repr(float(anchor @ vectors[relation["positive_domain_id"]])),
            "negative_score": repr(float(anchor @ vectors[relation["negative_domain_id"]])),
        })
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
