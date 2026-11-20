from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class LadderScore:
    ladder_id: str
    first_edit_event_id: str
    parent_locus_id: str
    outer_cluster_id: str
    cosine_x0_x1: float
    cosine_x0_x2: float
    cosine_x2_x1: float
    cosine_x2_x0: float
    delta_left: float
    delta_right: float
    F_credit_left: float
    F_credit_right: float
    ladder_F: float


def safe_cosine(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if left.ndim != 1 or right.ndim != 1 or left.shape != right.shape:
        raise ValueError("vectors must be one-dimensional with identical shape")
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("non-finite representation")
    left_scale = float(np.max(np.abs(left)))
    right_scale = float(np.max(np.abs(right)))
    if left_scale == 0.0 or right_scale == 0.0:
        return 0.0
    left_unit = left / left_scale
    right_unit = right / right_scale
    left_norm = float(np.linalg.norm(left_unit))
    right_norm = float(np.linalg.norm(right_unit))
    value = float(np.dot(left_unit / left_norm, right_unit / right_norm))
    return float(np.clip(value, -1.0, 1.0))


def contrast_F_credit(delta: float, tau: float) -> float:
    if not np.isfinite(delta) or not np.isfinite(tau) or tau < 0:
        raise ValueError("delta must be finite and tau must be finite and non-negative")
    if delta > tau:
        return 1.0
    if delta < -tau:
        return -1.0
    return 0.0


def score_ladder(
    row: Mapping[str, str], embeddings: Mapping[str, np.ndarray], tau: float
) -> LadderScore:
    x0 = embeddings[row["x0_object_id"]]
    x1 = embeddings[row["x1_object_id"]]
    x2 = embeddings[row["x2_object_id"]]
    s01 = safe_cosine(x0, x1)
    s02 = safe_cosine(x0, x2)
    s21 = safe_cosine(x2, x1)
    s20 = safe_cosine(x2, x0)
    delta_left = s01 - s02
    delta_right = s21 - s20
    credit_left = contrast_F_credit(delta_left, tau)
    credit_right = contrast_F_credit(delta_right, tau)
    return LadderScore(
        ladder_id=row["ladder_id"],
        first_edit_event_id=row["first_edit_event_id"],
        parent_locus_id=row["parent_locus_id"],
        outer_cluster_id=row["outer_cluster_id"],
        cosine_x0_x1=s01,
        cosine_x0_x2=s02,
        cosine_x2_x1=s21,
        cosine_x2_x0=s20,
        delta_left=delta_left,
        delta_right=delta_right,
        F_credit_left=credit_left,
        F_credit_right=credit_right,
        ladder_F=(credit_left + credit_right) / 2.0,
    )


def _mean_by(items: Iterable[tuple[str, float]]) -> dict[str, float]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for key, value in items:
        grouped[key].append(float(value))
    return {key: float(np.mean(values)) for key, values in grouped.items()}


def aggregate_hierarchy(
    scores: Sequence[LadderScore],
) -> tuple[float, dict[str, float], dict[str, float], dict[str, float]]:
    if not scores:
        raise ValueError("no ladder scores")
    query = _mean_by((s.first_edit_event_id, s.ladder_F) for s in scores)
    query_meta: dict[str, tuple[str, str]] = {}
    for score in scores:
        meta = (score.parent_locus_id, score.outer_cluster_id)
        previous = query_meta.setdefault(score.first_edit_event_id, meta)
        if previous != meta:
            raise ValueError("first-edit query maps to multiple loci or clusters")
    locus = _mean_by((query_meta[q][0], value) for q, value in query.items())
    locus_cluster: dict[str, str] = {}
    for locus_id, cluster_id in query_meta.values():
        previous = locus_cluster.setdefault(locus_id, cluster_id)
        if previous != cluster_id:
            raise ValueError("locus maps to multiple clusters")
    cluster = _mean_by((locus_cluster[locus_id], value) for locus_id, value in locus.items())
    return float(np.mean(list(cluster.values()))), query, locus, cluster


def bootstrap_cluster_distribution(
    cluster_scores: Mapping[str, float], draws: int, seed: int
) -> np.ndarray:
    if draws <= 0:
        raise ValueError("draws must be positive")
    keys = sorted(cluster_scores)
    if not keys:
        raise ValueError("no cluster scores")
    values = np.asarray([cluster_scores[key] for key in keys], dtype=np.float64)
    rng = np.random.default_rng(seed)
    distribution = np.empty(draws, dtype=np.float64)
    for draw in range(draws):
        indices = rng.integers(0, len(values), len(values))
        distribution[draw] = float(values[indices].mean())
    return distribution


def bootstrap_cluster_macro(
    cluster_scores: Mapping[str, float], draws: int, seed: int
) -> tuple[float, float]:
    distribution = bootstrap_cluster_distribution(cluster_scores, draws, seed)
    low, high = np.quantile(distribution, [0.025, 0.975])
    return float(low), float(high)


