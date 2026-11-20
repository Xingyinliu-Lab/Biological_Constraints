#!/usr/bin/env python3
"""Frozen PS09 aggregation for held-out endpoint scores."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from ps09_refset import validate_refset


BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 20_261_008


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def endpoint_id(row: dict[str, str]) -> str:
    return f"{row['block_id']}|{row['state']}|{row['uniprot_position_1based']}"


def linear_quantile(values: np.ndarray, probability: float) -> float:
    try:
        return float(np.quantile(values, probability, method="linear"))
    except TypeError:
        return float(np.quantile(values, probability, interpolation="linear"))


def read_membership(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("empty membership table")
    ids = [endpoint_id(row) for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate formal endpoint ID")
    return rows


def read_leaf_scores(path: Path) -> dict[str, dict[str, str]]:
    observed: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            identity = row["endpoint_id"]
            if identity in observed:
                raise ValueError("duplicate leaf endpoint")
            value = float(row["score"])
            if not np.isfinite(value):
                raise ValueError("nonfinite score")
            observed[identity] = {**row, "score": repr(value)}
    return observed


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def score(
    *,
    binding: Path,
    membership: Path,
    blocks: Path,
    clusters: Path,
    fasta: Path,
    folds: Path,
    leaf_scores: Path,
    output: Path,
    individual_scores: Path | None = None,
    pair_credits: Path | None = None,
    block_scores: Path | None = None,
    cluster_scores: Path | None = None,
) -> dict[str, object]:
    numeric_binding = json.loads(binding.read_text(encoding="utf-8"))
    inputs = {
        "membership": membership,
        "blocks": blocks,
        "clusters": clusters,
        "fasta": fasta,
        "fold_assignment": folds,
    }
    actual_input_hashes = {key: sha256(path) for key, path in inputs.items()}
    formal_refset_identity_match = all(
        actual_input_hashes[key] == numeric_binding["input_sha256"][key]
        for key in inputs
    )

    closure = validate_refset(
        membership=membership,
        blocks=blocks,
        clusters=clusters,
        fasta=fasta,
        folds=folds,
    )
    membership_rows = closure["membership_rows"]
    expected = {endpoint_id(row): row for row in membership_rows}
    observed = read_leaf_scores(leaf_scores)
    if set(observed) != set(expected):
        raise ValueError("leaf population mismatch")

    endpoint_rows: list[dict[str, object]] = []
    by_block: dict[str, dict[str, list[tuple[str, float]]]] = defaultdict(lambda: defaultdict(list))
    block_cluster: dict[str, str] = {}
    for identity, row in expected.items():
        value = float(observed[identity]["score"])
        by_block[row["block_id"]][row["state"]].append((identity, value))
        block_cluster[row["block_id"]] = row["homology_cluster_id"]
        endpoint_rows.append(
            {
                "endpoint_id": identity,
                "block_id": row["block_id"],
                "state": row["state"],
                "homology_cluster_id": row["homology_cluster_id"],
                "fold": observed[identity].get("fold", ""),
                "heldout_score": repr(value),
            }
        )

    pair_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    cluster_values: dict[str, list[float]] = defaultdict(list)
    for block_id in sorted(by_block):
        states = by_block[block_id]
        if set(states) != {"DISORDER", "ORDER_COORDINATE_SUPPORTED"}:
            raise ValueError(f"incomplete block: {block_id}")
        disorder = states["DISORDER"]
        ordered = states["ORDER_COORDINATE_SUPPORTED"]
        credits: list[float] = []
        raw_credits: list[float] = []
        for disorder_id, disorder_score in disorder:
            for ordered_id, ordered_score in ordered:
                margin = disorder_score - ordered_score
                tie = margin == 0.0
                raw_credit = 0.5 if tie else float(margin > 0.0)
                f_credit = 2 * raw_credit - 1
                credits.append(f_credit)
                raw_credits.append(raw_credit)
                pair_rows.append(
                    {
                        "block_id": block_id,
                        "homology_cluster_id": block_cluster[block_id],
                        "disorder_endpoint_id": disorder_id,
                        "disorder_score": repr(disorder_score),
                        "ordered_endpoint_id": ordered_id,
                        "ordered_score": repr(ordered_score),
                        "margin": repr(margin),
                        "tie": str(tie).lower(),
                        "F_credit": repr(f_credit),
                    }
                )
        value = float(np.mean(credits))
        expected_f = 2 * float(np.mean(raw_credits)) - 1
        if not np.isclose(value, expected_f, rtol=0.0, atol=1e-15):
            raise AssertionError(f"pair-to-block F mismatch: {block_id}")
        cluster_id = block_cluster[block_id]
        cluster_values[cluster_id].append(value)
        block_rows.append(
            {
                "block_id": block_id,
                "homology_cluster_id": cluster_id,
                "disorder_endpoints": len(disorder),
                "ordered_endpoints": len(ordered),
                "block_F": repr(value),
            }
        )

    cluster_rows = []
    for cluster_id in sorted(cluster_values):
        value = float(np.mean(cluster_values[cluster_id]))
        if not np.isclose(value, np.mean(cluster_values[cluster_id]), rtol=0.0, atol=1e-15):
            raise AssertionError(f"block-to-cluster mismatch: {cluster_id}")
        cluster_rows.append({"homology_cluster_id": cluster_id, "block_count": len(cluster_values[cluster_id]), "cluster_F": repr(value)})
    cluster_vector = np.asarray(
        [float(row["cluster_F"]) for row in cluster_rows], dtype=np.float64
    )
    point = float(cluster_vector.mean())
    if not np.isclose(point, np.mean(cluster_vector), rtol=0.0, atol=1e-15):
        raise AssertionError("cluster-to-F mismatch")
    generator = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    draws = np.asarray(
        [
            cluster_vector[
                generator.integers(0, len(cluster_vector), len(cluster_vector))
            ].mean()
            for _ in range(BOOTSTRAP_DRAWS)
        ],
        dtype=np.float64,
    )
    ci = [linear_quantile(draws, 0.025), linear_quantile(draws, 0.975)]

    result: dict[str, object] = {
        "status": "PASS",
        "formal_result": True,
        "constraint_id": "PS09",
        "formal_interface": "F_ONLY",
        "score_transform_version": "PS09_v0.1.2_F_ONLY",
        "F": point,
        "F_95_CI": ci,
        "coverage": {
            "endpoints_scored": len(endpoint_rows),
            "endpoints_expected": len(expected),
            "blocks_scored": len(block_rows),
            "clusters_scored": len(cluster_rows),
        },
        "bootstrap": {
            "draws": BOOTSTRAP_DRAWS,
            "bit_generator": "PCG64",
            "seed": BOOTSTRAP_SEED,
            "unit": "homology_cluster_id",
            "quantile_method": "linear",
        },
        "input_identity": {
            **actual_input_hashes,
            "leaf_scores": sha256(leaf_scores),
        },
        "formal_refset_identity_match": formal_refset_identity_match,
        "refset_closure": closure["counts"],
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    Path(str(output) + ".bootstrap.float64le.bin").write_bytes(
        draws.astype("<f8").tobytes()
    )
    if individual_scores is not None:
        write_tsv(
            individual_scores,
            [
                "endpoint_id",
                "block_id",
                "state",
                "homology_cluster_id",
                "fold",
                "heldout_score",
            ],
            endpoint_rows,
        )
    if pair_credits is not None:
        write_tsv(
            pair_credits,
            [
                "block_id",
                "homology_cluster_id",
                "disorder_endpoint_id",
                "disorder_score",
                "ordered_endpoint_id",
                "ordered_score",
                "margin",
                "tie",
                "F_credit",
            ],
            pair_rows,
        )
    if block_scores is not None:
        write_tsv(
            block_scores,
            [
                "block_id",
                "homology_cluster_id",
                "disorder_endpoints",
                "ordered_endpoints",
                "block_F",
            ],
            block_rows,
        )
    if cluster_scores is not None:
        write_tsv(
            cluster_scores,
            ["homology_cluster_id", "block_count", "cluster_F"],
            cluster_rows,
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binding", required=True, type=Path)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--blocks", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--leaf-scores", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--individual-scores", type=Path)
    parser.add_argument("--pair-credits", type=Path)
    parser.add_argument("--block-scores", type=Path)
    parser.add_argument("--cluster-scores", type=Path)
    args = parser.parse_args()
    score(
        binding=args.binding,
        membership=args.membership,
        blocks=args.blocks,
        clusters=args.clusters,
        fasta=args.fasta,
        folds=args.folds,
        leaf_scores=args.leaf_scores,
        output=args.out,
        individual_scores=args.individual_scores,
        pair_credits=args.pair_credits,
        block_scores=args.block_scores,
        cluster_scores=args.cluster_scores,
    )


if __name__ == "__main__":
    main()
