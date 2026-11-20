"""DNA03 public core-local tau recomputation."""
from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from pathlib import Path

import numpy as np

FLOOR = 1e-8
MAX_AUTHORIZED_TAU = 1e-4
TRANSITIONS = {
    "core100_to_250",
    "core100_to_500",
    "core100_to_1000",
    "core100_to_2000",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cosine_distance(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    denominator = np.linalg.norm(a) * np.linalg.norm(b)
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("invalid cosine vector")
    return float(1 - np.dot(a, b) / denominator)


def load_embeddings(path: Path):
    data = np.load(path, allow_pickle=False)
    if set(data.files) != {"view_ids", "embeddings"}:
        raise ValueError("embeddings NPZ must contain exactly view_ids and embeddings")
    ids = [str(value) for value in data["view_ids"]]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate view ID")
    raw = data["embeddings"]
    if (
        not np.issubdtype(raw.dtype, np.number)
        or np.issubdtype(raw.dtype, np.complexfloating)
        or np.issubdtype(raw.dtype, np.bool_)
    ):
        raise ValueError("embeddings must use a real numeric dtype")
    array = np.asarray(raw, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] <= 0 or array.shape[0] != len(ids) or not np.isfinite(array).all():
        raise ValueError("invalid core-local embeddings")
    return {value: index for index, value in enumerate(ids)}, array


def compute_tau(controls_path: Path, primary_path: Path, repeat_path: Path):
    primary_index, primary = load_embeddings(primary_path)
    repeat_index, repeat = load_embeddings(repeat_path)
    if primary_index.keys() != repeat_index.keys():
        raise ValueError("primary/repeat view universes differ")
    with controls_path.open(encoding="utf-8", newline="") as handle:
        controls = list(csv.DictReader(handle, delimiter="\t"))
    if len(controls) != 128 or len({row["query_id"] for row in controls}) != 128:
        raise ValueError("expected 128 unique controls")
    by_transition = defaultdict(list)
    ledger = []
    for row in controls:
        query, positive, negative = (
            row["query_view_id"],
            row["true_mate_view_id"],
            row["negative_view_id"],
        )
        if any(value not in primary_index for value in (query, positive, negative)):
            raise ValueError("control view missing")
        primary_d_negative = cosine_distance(primary[primary_index[query]], primary[primary_index[negative]])
        primary_d_positive = cosine_distance(primary[primary_index[query]], primary[primary_index[positive]])
        repeat_d_negative = cosine_distance(repeat[repeat_index[query]], repeat[repeat_index[negative]])
        repeat_d_positive = cosine_distance(repeat[repeat_index[query]], repeat[repeat_index[positive]])
        primary_contrast = primary_d_negative - primary_d_positive
        repeat_contrast = repeat_d_negative - repeat_d_positive
        delta = abs(primary_contrast - repeat_contrast)
        by_transition[row["transition_id"]].append(delta)
        ledger.append(
            {
                "query_id": row["query_id"],
                "transition_id": row["transition_id"],
                "primary_d_positive": primary_d_positive,
                "primary_d_negative": primary_d_negative,
                "primary_contrast": primary_contrast,
                "repeat_d_positive": repeat_d_positive,
                "repeat_d_negative": repeat_d_negative,
                "repeat_contrast": repeat_contrast,
                "abs_repeat_contrast_delta": delta,
            }
        )
    if set(by_transition) != TRANSITIONS or any(len(values) != 32 for values in by_transition.values()):
        raise ValueError("controls must contain 32 per transition")
    strata = {
        key: {
            "n": len(values),
            "max": float(np.max(values)),
            "q99": float(np.quantile(values, 0.99)),
        }
        for key, values in sorted(by_transition.items())
    }
    return {
        "selected_tau": max(FLOOR, max(item["q99"] for item in strata.values())),
        "maximum_abs_repeat_contrast_delta": max(item["max"] for item in strata.values()),
        "strata": strata,
        "ledger": ledger,
    }


def verify_tau_evidence(evidence: dict, controls_path: Path, primary_path: Path, repeat_path: Path):
    actual = compute_tau(controls_path, primary_path, repeat_path)
    if evidence.get("status") != "PASS" or evidence.get("controls") != 128:
        raise ValueError("invalid tau evidence status/count")
    expected = actual
    if abs(float(evidence["selected_tau"]) - expected["selected_tau"]) > 1e-15:
        raise ValueError("tau selected value mismatch")
    if evidence.get("ledger") != expected["ledger"]:
        raise ValueError("tau control ledger mismatch")
    if evidence.get("strata") != expected["strata"]:
        raise ValueError("tau strata mismatch")
    if expected["maximum_abs_repeat_contrast_delta"] > MAX_AUTHORIZED_TAU:
        raise ValueError("unstable repeat interface")
    return actual
