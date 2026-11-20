#!/usr/bin/env python3
"""Frozen RNA10 weighted logistic-ridge cross-fit measurement core."""

from __future__ import annotations

import hashlib
from collections import defaultdict

import numpy as np
from scipy.optimize import minimize


RIDGE_GRID = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0]


def mean(values):
    if not values:
        raise ValueError("empty value collection")
    return float(np.mean(values))


def weights(rows):
    groups = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for index, row in enumerate(rows):
        groups[row["source_dependency_unit_id"]][row["matched_block_id"]][row["truth_class"]].append(index)
    result = np.zeros(len(rows))
    for blocks in groups.values():
        for classes in blocks.values():
            for indices in classes.values():
                result[indices] = 1 / len(groups) / len(blocks) / len(classes) / len(indices)
    return result / result.mean()


def fit(matrix, rows, train, test, ridge_lambda):
    sample_weights = weights([rows[index] for index in train])
    sample_weights = sample_weights / sample_weights.sum()
    center = (matrix[train] * sample_weights[:, None]).sum(0)
    scale = np.sqrt(((matrix[train] - center) ** 2 * sample_weights[:, None]).sum(0))
    scale[scale == 0] = 1
    train_matrix = (matrix[train] - center) / scale
    test_matrix = (matrix[test] - center) / scale
    labels = np.array([1.0 if rows[index]["truth_class"] == "DISRUPTING" else 0.0 for index in train])

    def objective(parameters):
        logits = train_matrix @ parameters[:-1] + parameters[-1]
        probabilities = 1 / (1 + np.exp(-np.clip(logits, -40, 40)))
        residual = (probabilities - labels) * sample_weights
        loss = np.sum(sample_weights * (np.logaddexp(0, logits) - labels * logits))
        loss += 0.5 * ridge_lambda * np.dot(parameters[:-1], parameters[:-1])
        gradient = np.r_[
            train_matrix.T @ residual + ridge_lambda * parameters[:-1],
            residual.sum(),
        ]
        return float(loss), gradient

    fitted = minimize(
        lambda value: objective(value),
        np.zeros(train_matrix.shape[1] + 1),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": 2000},
    )
    if not fitted.success:
        raise RuntimeError(fitted.message)
    return test_matrix @ fitted.x[:-1] + fitted.x[-1]


def relation_score(rows, scores):
    blocks = defaultdict(list)
    for row, score in zip(rows, scores):
        blocks[row["matched_block_id"]].append((row, float(score)))
    source_scores = defaultdict(list)
    for values in blocks.values():
        disrupting = [score for row, score in values if row["truth_class"] == "DISRUPTING"]
        preserving = [score for row, score in values if row["truth_class"] == "PRESERVING"]
        if not disrupting or not preserving:
            raise ValueError("incomplete matched block")
        credit = mean([1.0 if left > right else 0.0 if left < right else 0.5 for left in disrupting for right in preserving])
        source_scores[values[0][0]["source_dependency_unit_id"]].append(credit)
    return mean([mean(values) for values in source_scores.values()])


def build_measurement_records(rows, matrix, plan):
    units = np.array([row["source_dependency_unit_id"] for row in rows])
    predictions = np.full(len(rows), np.nan)
    outer = []
    for held_out in sorted(set(units)):
        train = np.where(units != held_out)[0]
        test = np.where(units == held_out)[0]
        folds = defaultdict(list)
        for unit in sorted(set(units[train])):
            fold = int(hashlib.sha256(unit.encode()).hexdigest()[:8], 16) % 5
            folds[fold].append(unit)
        best = None
        for ridge_lambda in RIDGE_GRID:
            scores = []
            for fold_units in folds.values():
                validation = train[np.isin(units[train], fold_units)]
                inner_train = train[~np.isin(units[train], fold_units)]
                if len(validation) and len(inner_train):
                    scores.append(
                        relation_score(
                            [rows[index] for index in validation],
                            fit(matrix, rows, inner_train, validation, ridge_lambda),
                        )
                    )
            key = (mean(scores), -ridge_lambda)
            if best is None or key > best[0]:
                best = (key, ridge_lambda)
        predictions[test] = fit(matrix, rows, train, test, best[1])
        outer.append(
            {
                "held_out_source_dependency_unit_id": held_out,
                "selected_lambda": best[1],
                "inner_R": best[0][0],
            }
        )
    if not np.isfinite(predictions).all():
        raise ValueError("incomplete held-out prediction")

    blocks = defaultdict(lambda: defaultdict(list))
    for row, prediction in zip(rows, predictions):
        blocks[row["matched_block_id"]][row["truth_class"]].append((row, float(prediction)))
    records = []
    for block, arms in blocks.items():
        for disrupting, h_positive in arms["DISRUPTING"]:
            for preserving, h_negative in arms["PRESERVING"]:
                records.append(
                    {
                        "leaf_id": f"{block}::D:{disrupting['relation_id']}::P:{preserving['relation_id']}",
                        "matched_block_id": block,
                        "source_dependency_unit_id": disrupting["source_dependency_unit_id"],
                        "uncertainty_dependency_cluster_id": disrupting["uncertainty_dependency_cluster_id"],
                        "h_positive": h_positive,
                        "h_negative": h_negative,
                        "source_relation_ids": {
                            "disrupting": disrupting["relation_id"],
                            "preserving": preserving["relation_id"],
                        },
                    }
                )
    expected = {member["leaf_id"] for member in plan["canonical_members"]}
    if {record["leaf_id"] for record in records} != expected:
        raise ValueError("measurement membership mismatch")
    return sorted(records, key=lambda record: record["leaf_id"]), outer
