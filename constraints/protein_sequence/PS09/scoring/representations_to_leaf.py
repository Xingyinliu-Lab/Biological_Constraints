"""Convert third-party PS09 endpoint representations into frozen cross-fit leaf scores."""

from __future__ import annotations

import argparse
import csv
import gzip
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression

from ps09_refset import validate_refset


def endpoint_id(row: dict[str, str]) -> str:
    return f"{row['block_id']}|{row['state']}|{row['uniprot_position_1based']}"


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


def read_folds(path: Path) -> dict[str, int]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    folds = {row["homology_cluster_id"]: int(row["fold"]) for row in rows}
    if len(folds) != len(rows) or set(folds.values()) != set(range(5)):
        raise ValueError("invalid frozen five-fold assignment")
    return folds


def load_representations(path: Path, expected_ids: set[str]) -> dict[str, np.ndarray]:
    if path.suffix == ".safetensors":
        try:
            from safetensors import safe_open
        except ImportError as exc:
            raise RuntimeError(
                "safetensors input requires the safetensors package"
            ) from exc
        with safe_open(str(path), framework="np") as handle:
            ids = list(handle.keys())
            if set(ids) != expected_ids:
                missing = sorted(expected_ids - set(ids))[:5]
                extra = sorted(set(ids) - expected_ids)[:5]
                raise ValueError(
                    f"representation ID closure failure; missing={missing}, extra={extra}"
                )
            vectors = {}
            for identity in ids:
                raw = np.asarray(handle.get_tensor(identity))
                if not np.issubdtype(raw.dtype, np.number) or np.issubdtype(
                    raw.dtype, np.complexfloating
                ):
                    raise ValueError("representations must use a real numeric dtype")
                vectors[identity] = np.asarray(raw, dtype=np.float64)
        dimensions = {vector.shape for vector in vectors.values()}
        if len(dimensions) != 1 or not dimensions or len(next(iter(dimensions))) != 1:
            raise ValueError("each safetensors value must be one fixed-length vector")
        if next(iter(dimensions))[0] == 0:
            raise ValueError("embedding dimension must be positive")
        if not all(np.isfinite(vector).all() for vector in vectors.values()):
            raise ValueError("representations contain NaN or Inf")
        return vectors
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {"object_ids", "vectors"}:
            raise ValueError("NPZ must contain exactly object_ids and vectors")
        ids = [str(value) for value in data["object_ids"].tolist()]
        raw_vectors = np.asarray(data["vectors"])
        if not np.issubdtype(raw_vectors.dtype, np.number) or np.issubdtype(
            raw_vectors.dtype, np.complexfloating
        ):
            raise ValueError("representations must use a real numeric dtype")
        vectors = np.asarray(raw_vectors, dtype=np.float64)
    if vectors.ndim != 2 or vectors.shape[0] != len(ids) or vectors.shape[1] == 0:
        raise ValueError("vectors must have shape [n_objects, embedding_dim]")
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate object ID in representations")
    if set(ids) != expected_ids:
        missing = sorted(expected_ids - set(ids))[:5]
        extra = sorted(set(ids) - expected_ids)[:5]
        raise ValueError(f"representation ID closure failure; missing={missing}, extra={extra}")
    if not np.isfinite(vectors).all():
        raise ValueError("representations contain NaN or Inf")
    return dict(zip(ids, vectors, strict=True))


def validate_training_weight_closure(
    rows: list[dict[str, str]], train_indices: np.ndarray, weights: np.ndarray, fold: int
) -> float:
    training_states_by_block: dict[str, set[str]] = defaultdict(set)
    for index in train_indices:
        training_states_by_block[rows[index]["block_id"]].add(rows[index]["state"])
    if any(
        states != {"ORDER_COORDINATE_SUPPORTED", "DISORDER"}
        for states in training_states_by_block.values()
    ):
        raise ValueError(f"incomplete block-state training closure for fold {fold}")
    expected_total_weight = float(len(training_states_by_block))
    total_weight = float(weights[train_indices].sum())
    if not np.isclose(total_weight, expected_total_weight, rtol=0.0, atol=1e-12):
        raise ValueError(
            f"block-state training-weight closure failed for fold {fold}: "
            f"observed={total_weight}, expected={expected_total_weight}"
        )
    return expected_total_weight


def crossfit(
    rows: list[dict[str, str]], folds: dict[str, int], vectors: dict[str, np.ndarray]
) -> np.ndarray:
    if {row["homology_cluster_id"] for row in rows} != set(folds):
        raise ValueError("membership/fold cluster closure failure")
    matrix = np.stack([vectors[endpoint_id(row)] for row in rows])
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        groups[(row["block_id"], row["state"])].append(index)

    scores = np.full(len(rows), np.nan, dtype=np.float64)
    for fold in range(5):
        held_out = np.array(
            [folds[row["homology_cluster_id"]] == fold for row in rows], dtype=bool
        )
        train = ~held_out
        if not held_out.any() or not train.any():
            raise ValueError(f"empty train or held-out population for fold {fold}")

        weights = np.zeros(len(rows), dtype=np.float64)
        for indices in groups.values():
            selected = np.asarray(indices, dtype=int)
            selected = selected[train[selected]]
            if selected.size:
                weights[selected] = 1.0 / (2.0 * selected.size)
        train_indices = np.where(train)[0]
        test_indices = np.where(held_out)[0]
        if not np.all(weights[train_indices] > 0):
            raise ValueError("invalid block-state training weights")

        total_weight = validate_training_weight_closure(rows, train_indices, weights, fold)
        mean = (
            matrix[train_indices] * weights[train_indices, None]
        ).sum(axis=0) / total_weight
        variance = (
            weights[train_indices, None] * (matrix[train_indices] - mean) ** 2
        ).sum(axis=0) / total_weight
        scale = np.sqrt(variance)
        scale[scale == 0] = 1.0
        labels = np.asarray(
            [1 if rows[index]["state"] == "DISORDER" else 0 for index in train_indices],
            dtype=int,
        )
        if set(labels.tolist()) != {0, 1}:
            raise ValueError(f"fold {fold} training labels do not contain both classes")
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            classifier = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000)
            classifier.fit(
                (matrix[train_indices] - mean) / scale,
                labels,
                sample_weight=weights[train_indices],
            )
        if int(np.max(classifier.n_iter_)) >= 1000:
            raise RuntimeError(f"logistic regression did not converge for fold {fold}")
        scores[test_indices] = classifier.predict_proba(
            (matrix[test_indices] - mean) / scale
        )[:, 1]
    if not np.isfinite(scores).all():
        raise RuntimeError("cross-fit produced missing or non-finite scores")
    return scores


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--blocks", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--representations", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    closure = validate_refset(
        membership=args.membership,
        blocks=args.blocks,
        clusters=args.clusters,
        fasta=args.fasta,
        folds=args.folds,
    )
    rows = closure["membership_rows"]
    folds = closure["folds"]
    expected_ids = {endpoint_id(row) for row in rows}
    vectors = load_representations(args.representations, expected_ids)
    scores = crossfit(rows, folds, vectors)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["endpoint_id", "score", "fold"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(
            {
                "endpoint_id": endpoint_id(row),
                "score": repr(float(scores[index])),
                "fold": folds[row["homology_cluster_id"]],
            }
            for index, row in enumerate(rows)
        )


if __name__ == "__main__":
    main()
