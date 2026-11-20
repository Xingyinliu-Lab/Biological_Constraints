"""Frozen PS01 hierarchy and clan-bootstrap implementation."""
from __future__ import annotations

import csv
import gzip
import hashlib
import random
from collections import defaultdict
from pathlib import Path


LEVELS = (
    ("clan_id", "anchor_family", "anchor", "positive", "negative_family"),
    ("clan_id", "anchor_family", "anchor", "positive"),
    ("clan_id", "anchor_family", "anchor"),
    ("clan_id", "anchor_family"),
    ("clan_id",),
)


def average(values: list[float]) -> float:
    if not values:
        raise ValueError("cannot average an empty PS01 group")
    return sum(values) / len(values)


def read_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("PS01 contrast-credit ledger is empty")
    return rows


def hierarchy(rows: list[dict[str, str]]) -> tuple[float, dict[str, float]]:
    values = [({key: row[key] for key in LEVELS[0]}, float(row["F_credit"])) for row in rows]
    for keys in LEVELS:
        grouped: dict[tuple[str, ...], list[float]] = defaultdict(list)
        for metadata, value in values:
            grouped[tuple(metadata[key] for key in keys)].append(value)
        values = [
            ({key: value_key[index] for index, key in enumerate(keys)}, average(group_values))
            for value_key, group_values in grouped.items()
        ]
    clan_F = {metadata["clan_id"]: value for metadata, value in values}
    return average(list(clan_F.values())), clan_F


def evaluate(
    ledger: Path,
    bootstrap_replicates: int = 20_000,
    bootstrap_seed: int = 20_261_005,
    return_bootstrap: bool = False,
) -> dict:
    if bootstrap_replicates <= 0:
        raise ValueError("bootstrap_replicates must be positive")
    rows = read_rows(ledger)
    F, clan_F = hierarchy(rows)
    clans = sorted(clan_F)
    rng = random.Random(bootstrap_seed)
    F_bootstrap = [
        average([clan_F[rng.choice(clans)] for _ in clans])
        for _ in range(bootstrap_replicates)
    ]
    sorted_bootstrap = sorted(F_bootstrap)
    lower = sorted_bootstrap[int(0.025 * bootstrap_replicates)]
    upper = sorted_bootstrap[int(0.975 * bootstrap_replicates) - 1]
    F_credits = [float(row["F_credit"]) for row in rows]
    result = {
        "artifact_type": "PS01_F_ONLY_FORMAL_EVALUATOR_RESULT",
        "formal_interface": "F_ONLY",
        "score_transform_version": "PS01_v2.0.2_F_ONLY",
        "aggregation": "negative_family -> positive -> anchor -> anchor_family -> clan -> overall",
        "bootstrap": {
            "unit": "clan_id",
            "replicates": bootstrap_replicates,
            "seed": bootstrap_seed,
            "with_replacement": True,
            "draw_multiplicity_retained": True,
        },
        "row_count": len(rows),
        "leaf_wins": sum(value == 1.0 for value in F_credits),
        "leaf_ties": sum(value == 0.0 for value in F_credits),
        "leaf_losses": sum(value == -1.0 for value in F_credits),
        "leaf_mean_F": average(F_credits),
        "F": F,
        "F_95_CI": [lower, upper],
        "clan_F": clan_F,
        "ledger_sha256": hashlib.sha256(ledger.read_bytes()).hexdigest(),
    }
    if return_bootstrap:
        result["F_bootstrap"] = F_bootstrap
    return result
