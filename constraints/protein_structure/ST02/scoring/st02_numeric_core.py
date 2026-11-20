from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping

import numpy as np


def mean(values: Iterable[float]) -> float:
    values = list(values)
    if not values:
        raise ValueError("empty aggregation node")
    return float(np.mean(np.asarray(values, dtype=np.float64)))


def aggregate_hierarchy(
    rows: list[dict[str, str]],
    credits: list[float],
    cluster_weights: Mapping[str, float] | None = None,
) -> dict:
    if len(rows) != len(credits):
        raise ValueError("row/credit length mismatch")

    endpoint_values: dict[str, list[float]] = defaultdict(list)
    endpoint_meta: dict[str, tuple[str, str, str, str, str]] = {}
    for row, credit in zip(rows, credits):
        direction = row["direction_id"]
        endpoint_values[direction].append(float(credit))
        signature = (
            row["anchor_object_id"],
            row["target_homologous_superfamily_id"],
            row["source_dependency_cluster_id"],
            row["source_homologous_superfamily_id"],
            row["topology_id"],
        )
        if direction in endpoint_meta and endpoint_meta[direction] != signature:
            raise ValueError("direction metadata conflict")
        endpoint_meta[direction] = signature
    endpoint = {key: mean(values) for key, values in endpoint_values.items()}

    target_sf_values: dict[tuple, list[float]] = defaultdict(list)
    for direction, value in endpoint.items():
        target_sf_values[endpoint_meta[direction]].append(value)
    target_sf = {key: mean(values) for key, values in target_sf_values.items()}

    anchor_values: dict[tuple, list[float]] = defaultdict(list)
    for (anchor, target_sf_id, cluster, source_sf, topology), value in target_sf.items():
        anchor_values[(anchor, cluster, source_sf, topology)].append(value)
    anchor = {key: mean(values) for key, values in anchor_values.items()}

    cluster_values: dict[tuple, list[float]] = defaultdict(list)
    for (anchor_id, cluster, source_sf, topology), value in anchor.items():
        cluster_values[(cluster, source_sf, topology)].append(value)
    cluster = {key: mean(values) for key, values in cluster_values.items()}

    strata = sorted({(source_sf, topology) for _, source_sf, topology in cluster})
    source_sf_topology: dict[tuple[str, str], float] = {}
    for source_sf, topology in strata:
        values, weights = [], []
        for (cluster_id, sf_id, topology_id), value in cluster.items():
            if (sf_id, topology_id) != (source_sf, topology):
                continue
            weight = 1.0 if cluster_weights is None else float(cluster_weights[cluster_id])
            if not np.isfinite(weight) or weight <= 0:
                raise ValueError("cluster bootstrap weight must be finite and positive")
            values.append(value)
            weights.append(weight)
        source_sf_topology[(source_sf, topology)] = float(
            np.average(np.asarray(values, dtype=np.float64), weights=np.asarray(weights, dtype=np.float64))
        )

    topology: dict[str, float] = {}
    for topology_id in sorted({topology_id for _, topology_id in strata}):
        topology[topology_id] = mean(
            value for (source_sf, top), value in source_sf_topology.items() if top == topology_id
        )
    overall = mean(topology.values())
    return {
        "F": overall,
        "source_sf_topology": source_sf_topology,
        "topology": topology,
        "counts": {
            "endpoint": len(endpoint),
            "target_sf_within_anchor": len(target_sf),
            "anchor": len(anchor),
            "dependency_cluster_source_sf_topology": len(cluster),
            "source_sf_topology": len(source_sf_topology),
            "topology": len(topology),
        },
    }


def bayesian_cluster_bootstrap(
    rows: list[dict[str, str]],
    credits: list[float],
    *,
    draws: int,
    seed: int,
    return_overall_draws: bool = False,
) -> dict:
    clusters = sorted({row["source_dependency_cluster_id"] for row in rows})
    if not clusters:
        raise ValueError("no dependency clusters")
    rng = np.random.Generator(np.random.PCG64(seed))
    overall = np.empty(draws, dtype=np.float64)
    topology_values: dict[str, np.ndarray] = {}
    stratum_values: dict[tuple[str, str], np.ndarray] = {}
    for draw_index in range(draws):
        weights = rng.exponential(scale=1.0, size=len(clusters))
        result = aggregate_hierarchy(rows, credits, dict(zip(clusters, weights.tolist())))
        overall[draw_index] = result["F"]
        if not topology_values:
            topology_values = {key: np.empty(draws, dtype=np.float64) for key in result["topology"]}
            stratum_values = {key: np.empty(draws, dtype=np.float64) for key in result["source_sf_topology"]}
        for key, value in result["topology"].items():
            topology_values[key][draw_index] = value
        for key, value in result["source_sf_topology"].items():
            stratum_values[key][draw_index] = value

    def ci(values: np.ndarray) -> list[float]:
        return np.quantile(values, [0.025, 0.975], method="linear").astype(float).tolist()

    result = {
        "overall_ci95": ci(overall),
        "topology_ci95": {key: ci(values) for key, values in topology_values.items()},
        "source_sf_topology_ci95": {key: ci(values) for key, values in stratum_values.items()},
        "valid_draws": draws,
    }
    if return_overall_draws:
        result["overall_draws"] = overall
    return result
