"""Frozen PS10 cross-fit calculation shared by all public runs."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from ps10_refset import validate_refset


BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 20_261_008


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_fold_f(fold_f: list[float]) -> tuple[dict[str, object], np.ndarray]:
    fold_vector = np.asarray(fold_f, dtype=np.float64)
    if len(fold_vector) != 22 or not np.isfinite(fold_vector).all():
        raise ValueError("PS10 fold AUC closure failure")
    generator = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    draws = np.asarray(
        [
            generator.choice(fold_vector, len(fold_vector), replace=True).mean()
            for _ in range(BOOTSTRAP_DRAWS)
        ],
        dtype=np.float64,
    )
    return {
        "F": float(fold_vector.mean()),
        "F_95_CI": [
            float(np.quantile(draws, 0.025, method="linear")),
            float(np.quantile(draws, 0.975, method="linear")),
        ],
    }, draws


def evaluate(
    records_path: Path, ledger_path: Path, embeddings_path: Path
) -> dict[str, object]:
    closure = validate_refset(records_path, ledger_path)
    records = closure["records"]
    ledger = closure["assignments"]
    ids = [row["record_id"] for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate frozen record ID")
    by_id = {row["record_id"]: row for row in records}
    if len(ledger) != len(ids) * 22 or len(
        {(row["fold_id"], row["record_id"]) for row in ledger}
    ) != len(ledger):
        raise ValueError("frozen fold ledger invariant failed")

    with np.load(embeddings_path, allow_pickle=False) as archive:
        if set(archive.files) != {"record_id", "embedding"}:
            raise ValueError("NPZ must contain exactly record_id and embedding")
        embedding_ids = [str(value) for value in archive["record_id"]]
        # Preserve the executed representation dtype. The frozen paper runs
        # fitted the readout directly on the stored float32 embeddings; an
        # eager float64 cast changes a small number of rank ties and therefore
        # the held-out AUCs.
        matrix = np.asarray(archive["embedding"])
    if not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(
        matrix.dtype, np.complexfloating
    ):
        raise ValueError("embeddings must use a real numeric dtype")
    if matrix.ndim != 2 or matrix.shape[0] != len(embedding_ids) or matrix.shape[1] == 0:
        raise ValueError("invalid embedding matrix shape")
    if len(embedding_ids) != len(set(embedding_ids)) or set(embedding_ids) != set(ids):
        raise ValueError("embedding population mismatch")
    if not np.isfinite(matrix).all():
        raise ValueError("embeddings contain NaN or Inf")

    vectors = dict(zip(embedding_ids, matrix, strict=True))
    by_fold: dict[str, list[dict[str, str]]] = {}
    for row in ledger:
        if row["record_id"] not in by_id or row["role"] not in {
            "train",
            "test",
            "purged",
        }:
            raise ValueError("invalid frozen fold ledger row")
        by_fold.setdefault(row["fold_id"], []).append(row)

    fold_rows: list[dict[str, object]] = []
    individual_rows: list[dict[str, object]] = []
    fold_f: list[float] = []
    for fold_id, rows in sorted(by_fold.items()):
        train = [row["record_id"] for row in rows if row["role"] == "train"]
        test = [row["record_id"] for row in rows if row["role"] == "test"]
        purged = [row["record_id"] for row in rows if row["role"] == "purged"]
        if not train or not test or set(train) & set(test):
            raise ValueError(f"invalid train/test population in fold {fold_id}")

        studies = {
            study
            for record_id in train
            for study in by_id[record_id]["source_study_units"].split(";")
            if study
        }
        counts = {
            study: sum(
                record_id in train
                for record_id in ids
                if study in by_id[record_id]["source_study_units"].split(";")
            )
            for study in studies
        }
        weights = np.asarray(
            [
                sum(
                    1.0 / counts[study]
                    for study in by_id[record_id]["source_study_units"].split(";")
                    if study in counts
                )
                for record_id in train
            ],
            dtype=np.float64,
        )
        labels = np.asarray(
            [
                by_id[record_id]["experimental_classification"] == "AMYLOID"
                for record_id in train
            ],
            dtype=int,
        )
        training = np.stack([vectors[record_id] for record_id in train])
        mean = training.mean(axis=0)
        scale = training.std(axis=0)
        if np.any(scale == 0):
            raise ValueError("zero-SD training feature")

        model = LogisticRegression(
            C=1.0,
            solver="lbfgs",
            max_iter=1000,
            fit_intercept=True,
            class_weight=None,
        ).fit((training - mean) / scale, labels, sample_weight=weights)
        if int(model.n_iter_.max()) >= 1000:
            raise RuntimeError("nonconvergent PS10 fold fit")

        test_labels = np.asarray(
            [
                by_id[record_id]["experimental_classification"] == "AMYLOID"
                for record_id in test
            ],
            dtype=int,
        )
        decision_scores = model.decision_function(
            (np.stack([vectors[record_id] for record_id in test]) - mean) / scale
        )
        auc = float(roc_auc_score(test_labels, decision_scores))
        fold_score = 2.0 * auc - 1.0
        fold_f.append(fold_score)
        fold_rows.append(
            {
                "fold_id": fold_id,
                "train_records": len(train),
                "test_records": len(test),
                "purged_records": len(purged),
                "heldout_auc": repr(auc),
                "fold_F": repr(fold_score),
            }
        )
        for record_id, label, decision_score in zip(
            test, test_labels, decision_scores, strict=True
        ):
            individual_rows.append(
                {
                    "fold_id": fold_id,
                    "record_id": record_id,
                    "experimental_classification": by_id[record_id][
                        "experimental_classification"
                    ],
                    "binary_label": int(label),
                    "source_study_units": by_id[record_id]["source_study_units"],
                    "heldout_decision_score": repr(float(decision_score)),
                }
            )

    summary, draws = summarize_fold_f(fold_f)
    return {
        "summary": {
            **summary,
            "coverage": {
                "records": len(records),
                "folds": len(fold_rows),
                "heldout_prediction_rows": len(individual_rows),
            },
            "refset_closure": closure["counts"],
        },
        "individual_rows": individual_rows,
        "fold_rows": fold_rows,
        "bootstrap_draws": draws,
    }
