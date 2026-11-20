#!/usr/bin/env python3
"""Run the frozen DNA01 exact 4-mer reverse-complement reference."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

SCORING = Path(__file__).resolve().parents[3] / "scoring"
sys.path.insert(0, str(SCORING))
from dna01_core import (  # noqa: E402
    bounded_credit,
    publish_result,
    read_tsv,
    sha256_file,
    summarize_leaves,
    validate_refset,
)

ALPHABET = "ACGT"
KMERS = tuple(a + b + c + d for a in ALPHABET for b in ALPHABET for c in ALPHABET for d in ALPHABET)
KMER_INDEX = {kmer: index for index, kmer in enumerate(KMERS)}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
RC_PERMUTATION = np.asarray(
    [KMER_INDEX[kmer.translate(COMPLEMENT)[::-1]] for kmer in KMERS], dtype=np.int64
)


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise ValueError("duplicate FASTA identifier")
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence before identifier")
                chunks.append(line)
    if name is not None:
        if name in records:
            raise ValueError("duplicate FASTA identifier")
        records[name] = "".join(chunks).upper()
    return records


def fourmer_frequency(sequence):
    if len(sequence) < 4 or set(sequence) - set(ALPHABET):
        raise ValueError("sequence must contain only ACGT and have length at least four")
    vector = np.zeros(256, dtype=np.float64)
    for index in range(len(sequence) - 3):
        vector[KMER_INDEX[sequence[index:index + 4]]] += 1.0
    return vector / float(len(sequence) - 3)


def apply_rc(vector):
    output = np.empty_like(vector)
    output[RC_PERMUTATION] = vector
    return output


def main():
    parser = argparse.ArgumentParser(description="Run the DNA01 exact 4-mer RC white-box reference.")
    parser.add_argument("--relations", required=True, type=Path)
    parser.add_argument("--anchors", required=True, type=Path)
    parser.add_argument("--objects", required=True, type=Path)
    parser.add_argument("--fold-assignment", required=True, type=Path)
    parser.add_argument("--dependency-clusters", required=True, type=Path)
    parser.add_argument("--repeat-control-panel", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    relations = read_tsv(args.relations)
    anchors = read_tsv(args.anchors)
    objects = read_tsv(args.objects)
    folds = read_tsv(args.fold_assignment)
    dependency_rows = read_tsv(args.dependency_clusters)
    repeat_rows = read_tsv(args.repeat_control_panel)
    object_by_id, _relation_by_id, anchors_by_relation, _fold_by_cluster = validate_refset(
        objects, relations, anchors, folds, dependency_rows, repeat_rows
    )
    records = read_fasta(args.fasta)
    if not set(object_by_id).issubset(records):
        raise ValueError("FASTA is missing one or more formal objects")
    embeddings = {}
    for object_id, row in object_by_id.items():
        sequence = records[object_id]
        if len(sequence) != int(row["length_bp"]):
            raise ValueError(f"sequence length mismatch: {object_id}")
        if hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
            raise ValueError(f"sequence hash mismatch: {object_id}")
        embeddings[object_id] = fourmer_frequency(sequence)

    leaves = []
    for relation in relations:
        forward = embeddings[relation["forward_object_id"]]
        target = embeddings[relation["rc_object_id"]]
        predicted = apply_rc(forward)
        selected = anchors_by_relation[relation["relation_id"]]
        residual = float(np.linalg.norm(predicted - target))
        target_distances = [float(np.linalg.norm(target - embeddings[row["anchor_target_object_id"]])) for row in selected]
        target_scale = float(np.median(target_distances))
        wrong = [float(np.linalg.norm(predicted - embeddings[row["anchor_target_object_id"]])) for row in selected]
        cycle_residual = float(np.linalg.norm(apply_rc(predicted) - forward))
        cycle_scale = float(np.median([
            np.linalg.norm(forward - embeddings[row["anchor_forward_object_id"]]) for row in selected
        ]))
        credit = bounded_credit(residual, target_scale)
        leaf = {
            "relation_id": relation["relation_id"],
            "forward_object_id": relation["forward_object_id"],
            "rc_object_id": relation["rc_object_id"],
            "outer_fold_id": int(relation["outer_fold_id"]),
            "cluster_id": relation["source_dependency_cluster_id"],
            "parent_locus_id": relation["parent_locus_id"],
            "view_id": relation["view_id"],
            "evaluation_status": "SCORED" if credit is not None else "UNSCORABLE",
            "unscorable_reason": "" if credit is not None else "ZERO_TARGET_SCALE",
            "residual": residual,
            "target_scale": target_scale,
            "F_credit": credit,
            "specificity": float(np.mean([1.0 if residual < value else 0.5 if residual == value else 0.0 for value in wrong])),
            "cycle_residual": cycle_residual,
            "cycle_scale": cycle_scale,
            "normalized_cycle": cycle_residual / cycle_scale if cycle_scale > 0 else None,
        }
        for rank, distance in enumerate(wrong, start=1):
            leaf[f"anchor_target_distance_{rank}"] = distance
        for rank, distance in enumerate(target_distances, start=1):
            leaf[f"target_to_anchor_distance_{rank}"] = distance
        leaves.append(leaf)

    summary = summarize_leaves(leaves)
    result = {
        "status": "PASS",
        "score_status": summary["score_status"],
        "constraint_id": "C-DNA-01",
        "model_key": "DNA01_WB_4MER_RC_PERMUTATION",
        "result_kind": "ANALYTIC_WHITE_BOX_REFERENCE",
        "operator_branch": "EXACT_RC_4MER_CHANNEL_PERMUTATION",
        "score_transform_version": "DNA01_v0.5.4_F_ONLY_FULL_CHAIN",
        "score_transform": "F_credit = a_i/(a_i+e_i)",
        "F": summary["F"],
        "F_95_CI": summary["F_95_CI"],
        "Q_spec": summary["Q_spec"],
        "representation_status": "EXPLICIT_EQUIVARIANT_VALID",
        "equivariance_validity": "VALID",
        "coverage": summary["coverage"],
        "counts": {"objects": len(objects), "relations": len(relations), "anchors": len(anchors), "clusters": len(folds)},
        "operator_involution_exact": bool(np.array_equal(RC_PERMUTATION[RC_PERMUTATION], np.arange(256))),
        "input_sha256": {
            "relations": sha256_file(args.relations), "anchors": sha256_file(args.anchors),
            "objects": sha256_file(args.objects), "fold_assignment": sha256_file(args.fold_assignment),
            "dependency_clusters": sha256_file(args.dependency_clusters),
            "repeat_control_panel": sha256_file(args.repeat_control_panel), "fasta": sha256_file(args.fasta),
        },
    }
    publish_result(args.output, leaves, result, summary["parent_F"], summary["cluster_F"], summary["F_bootstrap"])
    print(json.dumps({"F": result["F"], "F_95_CI": result["F_95_CI"]}, sort_keys=True))


if __name__ == "__main__":
    main()
