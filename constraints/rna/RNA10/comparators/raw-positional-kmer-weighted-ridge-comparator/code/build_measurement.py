#!/usr/bin/env python3
"""Build the frozen RNA10 raw-kmer comparator measurement ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scoring"))
from rna10_measurement_core import build_measurement_records


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--effects", required=True, type=Path)
    parser.add_argument("--vectors", required=True, type=Path)
    parser.add_argument("--run-plan", required=True, type=Path)
    parser.add_argument("--binding", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--model-id", default="raw_positional_kmer_weighted_ridge")
    args = parser.parse_args()

    effects_sha256 = sha256(args.effects)
    vectors_sha256 = sha256(args.vectors)
    with args.effects.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    relation_ids = [row["relation_id"] for row in rows]
    if len(relation_ids) != len(set(relation_ids)):
        raise ValueError("duplicate relation ID")
    plan = json.loads(args.run_plan.read_text(encoding="utf-8"))
    binding = json.loads(args.binding.read_text(encoding="utf-8"))
    with np.load(args.vectors, allow_pickle=False) as archive:
        if set(archive.files) != {"sequence_sha256", "vectors"}:
            raise ValueError("vector NPZ must contain exactly sequence_sha256 and vectors")
        ids = [str(value) for value in archive["sequence_sha256"].tolist()]
        raw_vectors = archive["vectors"]
        if np.issubdtype(raw_vectors.dtype, np.bool_) or not np.issubdtype(raw_vectors.dtype, np.number) or np.issubdtype(raw_vectors.dtype, np.complexfloating):
            raise ValueError("vectors must use a real numeric dtype")
        vectors = np.asarray(raw_vectors, dtype=np.float64)
    if vectors.ndim != 2 or vectors.shape[0] != len(ids) or vectors.shape[1] <= 0:
        raise ValueError("invalid vector matrix shape")
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate sequence SHA256")
    vector_map = dict(zip(ids, vectors))
    required = {row["wt_sequence_sha256"] for row in rows} | {row["mut_sequence_sha256"] for row in rows}
    if set(vector_map) != required:
        raise ValueError("vector object universe mismatch")
    matrix = np.stack(
        [vector_map[row["mut_sequence_sha256"]] - vector_map[row["wt_sequence_sha256"]] for row in rows]
    ).astype(float)
    if not np.isfinite(matrix).all():
        raise ValueError("nonfinite representation delta")

    records, outer = build_measurement_records(rows, matrix, plan)
    output = {
        "artifact_type": "C_RNA10_MEASUREMENT_LEDGER",
        "model_id": args.model_id,
        "measurement_artifact_release_id": "RNA10_raw_positional_kmer_weighted_ridge_FORMAL_MEASUREMENT_ARTIFACT_v1.0.5_20261009",
        "runplan_sha256": hashlib.sha256(json.dumps(plan, sort_keys=True).encode()).hexdigest(),
        "population_identity_sha256": binding["population_identity_sha256"],
        "numeric_method_sha256": binding["numeric_method_sha256"],
        "split_release_sha256": binding["split_release_sha256"],
        "representation_dimension": int(matrix.shape[1]),
        "active_delta_feature_dimension": int(np.count_nonzero(np.any(matrix != 0, axis=0))),
        "feature_contract": "17_relative_position_bins_x_kmer_counts_k1_to_k6; runtime_score_uses_delta_phi_mut_minus_WT",
        "representation_semantics": "z is the raw positional k-mer feature vector phi; therefore z_mut-z_WT equals phi_mut-phi_WT",
        "delta_definition": "z_mut-z_WT",
        "label_encoding": {"PRESERVING": 0, "DISRUPTING": 1},
        "numeric_decoder": "weighted_logistic_ridge_lbfgsb",
        "source_effects_sha256": effects_sha256,
        "source_vectors_sha256": vectors_sha256,
        "source_inputs": {
            "vectors_sha256": vectors_sha256,
            "effects_sha256": effects_sha256,
            "runplan_file_sha256": sha256(args.run_plan),
            "binding_file_sha256": sha256(args.binding),
        },
        "records": records,
        "outer_LOO": outer,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
