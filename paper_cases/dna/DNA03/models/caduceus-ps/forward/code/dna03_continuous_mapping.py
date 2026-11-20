"""Frozen DNA03 v0.2.0 nucleotide-coordinate reconstruction primitives.

This module contains no model-specific code.  Model adapters must supply real
token spans (or native one-base spans) for each window.  Multi-base token spans
must never be guessed from token text.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


SUPPORTED_LENGTHS = (100, 250, 500, 1000, 2000)
CORE_LENGTH = 100


@dataclass(frozen=True)
class WindowProjection:
    start: int
    end: int
    positions: np.ndarray
    token_coverage: np.ndarray


def frozen_windows(length: int) -> tuple[tuple[int, int], ...]:
    """Return the only formal DNA03 v0.2.0 window layout."""
    if length not in SUPPORTED_LENGTHS:
        raise ValueError(f"unsupported formal DNA03 length: {length}")
    if length <= 1000:
        return ((0, length),)
    return ((0, 1000), (500, 1500), (1000, 2000))


def native_single_base_spans(length: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, i + 1) for i in range(length))


def project_tokens_to_nucleotides(
    hidden: np.ndarray,
    spans: Sequence[tuple[int, int]],
    nucleotide_length: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Average every valid token hidden state over each nucleotide it covers.

    Special tokens and padding must be represented by empty spans (0, 0) and
    are excluded.  All non-empty spans must be true tokenizer offsets within
    the current window.  Complete nucleotide coverage is mandatory.
    """
    h = np.asarray(hidden, dtype=np.float64)
    if h.ndim != 2 or h.shape[0] != len(spans):
        raise ValueError("hidden rows and token spans do not align")
    if not np.isfinite(h).all():
        raise ValueError("non-finite token hidden state")
    sums = np.zeros((nucleotide_length, h.shape[1]), dtype=np.float64)
    counts = np.zeros(nucleotide_length, dtype=np.int64)
    for row, (start, end) in zip(h, spans):
        start, end = int(start), int(end)
        if start == end == 0:
            continue
        if not (0 <= start < end <= nucleotide_length):
            raise ValueError(f"invalid token span [{start},{end})")
        sums[start:end] += row
        counts[start:end] += 1
    if np.any(counts == 0):
        missing = np.flatnonzero(counts == 0)
        raise ValueError(f"incomplete nucleotide coverage: {missing[:10].tolist()}")
    return (sums / counts[:, None]).astype(np.float32), counts


def merge_window_projections(
    sequence_length: int,
    projections: Iterable[WindowProjection],
) -> tuple[np.ndarray, np.ndarray]:
    """Map window-local positions to global coordinates and average overlaps."""
    projections = list(projections)
    observed = tuple((int(x.start), int(x.end)) for x in projections)
    expected = frozen_windows(sequence_length)
    if observed != expected:
        raise ValueError(f"window set/order mismatch: observed={observed}, expected={expected}")
    sums = None
    counts = np.zeros(sequence_length, dtype=np.int64)
    for projection in projections:
        start, end = int(projection.start), int(projection.end)
        arr = np.asarray(projection.positions, dtype=np.float64)
        if not (0 <= start < end <= sequence_length) or arr.shape[0] != end - start:
            raise ValueError("window projection coordinates or length are invalid")
        if np.any(np.asarray(projection.token_coverage) <= 0):
            raise ValueError("window has uncovered nucleotide positions")
        if sums is None:
            sums = np.zeros((sequence_length, arr.shape[1]), dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != sums.shape[1] or not np.isfinite(arr).all():
            raise ValueError("invalid window position array")
        sums[start:end] += arr
        counts[start:end] += 1
    if sums is None or np.any(counts == 0):
        raise ValueError("global nucleotide coverage is incomplete")
    return (sums / counts[:, None]).astype(np.float32), counts


def pool_core_local(
    position_states: np.ndarray,
    core_start: int,
    core_end: int,
) -> np.ndarray:
    """Return the frozen core-local nucleotide-coordinate mean."""
    arr = np.asarray(position_states, dtype=np.float64)
    if arr.ndim != 2 or not np.isfinite(arr).all():
        raise ValueError("invalid position-level representation")
    if core_end - core_start != CORE_LENGTH or not (0 <= core_start < core_end <= len(arr)):
        raise ValueError("formal core interval must contain exactly 100 nucleotides")
    return arr[core_start:core_end].mean(axis=0).astype(np.float32)


def require_real_offsets(
    offsets: Sequence[tuple[int, int]] | None,
    *,
    tokenizer_is_native_single_base: bool,
    nucleotide_length: int,
) -> tuple[tuple[int, int], ...]:
    """Fail closed when a multi-base tokenizer cannot return real offsets."""
    if offsets is not None:
        return tuple((int(a), int(b)) for a, b in offsets)
    if tokenizer_is_native_single_base:
        return native_single_base_spans(nucleotide_length)
    raise ValueError("multi-base tokenizer lacks authoritative token offsets")
