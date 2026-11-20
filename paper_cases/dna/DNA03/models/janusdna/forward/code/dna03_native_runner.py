"""Shared DNA03 v0.2.0 runner for native single-nucleotide hidden states."""
from __future__ import annotations

import csv, hashlib
from pathlib import Path

import numpy as np
import torch

from dna03_continuous_mapping import (
    WindowProjection, frozen_windows, merge_window_projections,
    native_single_base_spans, pool_core_local, project_tokens_to_nucleotides,
)


def read_fasta(path: Path):
    data, key = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]; data[key] = ""
        else: data[key] += line.strip().upper()
    return data


def read_rows(path: Path, sequences):
    with path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        seq = sequences[row["view_id"]]
        cs, ce = int(row["core_offset_start_0based"]), int(row["core_offset_end_exclusive"])
        if len(seq) != int(row["target_length_nt"]) or hashlib.sha256(seq.encode()).hexdigest() != row["view_sequence_sha256"] or ce - cs != 100:
            raise ValueError("view identity or core coordinates mismatch")
    return rows


def build_tasks(rows, sequences):
    tasks = []
    for row in rows:
        view, seq = row["view_id"], sequences[row["view_id"]]
        for start, end in frozen_windows(len(seq)):
            tasks.append((end-start, view, start, end, seq[start:end]))
    return sorted(tasks, key=lambda x: (x[0], x[1], x[2]))


def encode_native(rows, sequences, batch_size, forward_batch):
    tasks, by_view = build_tasks(rows, sequences), {}
    for base in range(0, len(tasks), batch_size):
        batch = tasks[base:base+batch_size]
        arrays = forward_batch([x[4] for x in batch])
        if len(arrays) != len(batch): raise ValueError("forward batch size mismatch")
        for (_, view, start, end, seq), hidden in zip(batch, arrays):
            hidden = np.asarray(hidden)
            if hidden.shape[0] != len(seq): raise ValueError("native hidden length mismatch")
            positions, coverage = project_tokens_to_nucleotides(hidden, native_single_base_spans(len(seq)), len(seq))
            by_view.setdefault(view, []).append(WindowProjection(start, end, positions, coverage))
    vectors, coverage_out = [], {}
    by_id = {r["view_id"]:r for r in rows}
    for view in [r["view_id"] for r in rows]:
        merged, counts = merge_window_projections(len(sequences[view]), by_view[view])
        row = by_id[view]
        vector = pool_core_local(merged, int(row["core_offset_start_0based"]), int(row["core_offset_end_exclusive"]))
        vectors.append(vector); coverage_out[view] = {"min":int(counts.min()),"max":int(counts.max())}
    return [r["view_id"] for r in rows], np.stack(vectors), len(tasks), coverage_out
