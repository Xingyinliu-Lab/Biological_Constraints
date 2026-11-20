"""Run the public frozen DNA02-02 scoring protocol on compatible representations."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from dna0202_numeric_core import (
    aggregate_hierarchy,
    bootstrap_cluster_distribution,
    bootstrap_cluster_macro,
    safe_cosine,
    score_ladder,
)
from dna0202_io import load_embeddings, read_tsv, validate_sequence_universe, write_tsv


DETERMINISTIC_FLOOR = 1e-8
MAX_AUTHORIZED_TAU = 1e-4
CONSTRAINT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFSET = CONSTRAINT_ROOT / "refset"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def contrast(row: dict[str, str], embeddings: dict[str, np.ndarray]) -> tuple[float, float]:
    x0, x1, x2 = (embeddings[row[key]] for key in ("x0_object_id", "x1_object_id", "x2_object_id"))
    return (
        safe_cosine(x0, x1) - safe_cosine(x0, x2),
        safe_cosine(x2, x1) - safe_cosine(x2, x0),
    )


def calibrate_tau(
    controls: list[dict[str, str]],
    ladders: list[dict[str, str]],
    primary: dict[str, np.ndarray],
    repeat: dict[str, np.ndarray],
) -> tuple[float, dict]:
    if set(primary) != set(repeat):
        raise ValueError("primary/repeat object universes differ")
    ladder_by_id = {row["ladder_id"]: row for row in ladders}
    if len(ladder_by_id) != len(ladders):
        raise ValueError("duplicate formal ladder ID")
    by_length: dict[int, list[float]] = defaultdict(list)
    seen: set[str] = set()
    for row in controls:
        if row["ladder_id"] in seen:
            raise ValueError("duplicate tau control ladder")
        seen.add(row["ladder_id"])
        formal = ladder_by_id.get(row["ladder_id"])
        if formal is None or any(row[key] != formal[key] for key in ("x0_object_id", "x1_object_id", "x2_object_id", "outer_cluster_id")):
            raise ValueError("tau control does not match its formal ladder")
        if any(row[key] not in primary for key in ("x0_object_id", "x1_object_id", "x2_object_id")):
            raise ValueError("tau control endpoint missing from embeddings")
        primary_values = contrast(row, primary)
        repeat_values = contrast(row, repeat)
        by_length[int(row["sequence_length"])].extend(
            abs(a - b) for a, b in zip(primary_values, repeat_values, strict=True)
        )
    expected_lengths = {100, 250, 500, 1000, 2000, 5000}
    if set(by_length) != expected_lengths or any(len(values) != 64 for values in by_length.values()):
        raise ValueError("tau controls must contain 32 ladders per frozen length stratum")
    strata = {}
    for length, values in sorted(by_length.items()):
        array = np.asarray(values, dtype=np.float64)
        strata[str(length)] = {
            "control_ladders": len(values) // 2,
            "endpoint_contrasts": len(values),
            "max_abs_repeat_delta": float(array.max()),
            "q99_abs_repeat_delta": float(np.quantile(array, 0.99)),
        }
    maximum = max(row["max_abs_repeat_delta"] for row in strata.values())
    tau = (
        DETERMINISTIC_FLOOR
        if maximum <= DETERMINISTIC_FLOOR
        else max(DETERMINISTIC_FLOOR, max(row["q99_abs_repeat_delta"] for row in strata.values()))
    )
    if not np.isfinite(tau) or tau > MAX_AUTHORIZED_TAU:
        raise ValueError(f"unstable repeat interface: selected tau={tau}")
    return tau, {
        "status": "PASS",
        "control_ladders": len(controls),
        "endpoint_contrasts": 2 * len(controls),
        "selected_tau": tau,
        "max_authorized_tau": MAX_AUTHORIZED_TAU,
        "maximum_abs_repeat_delta": maximum,
        "strata": strata,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Score model representations against DNA02-02.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ladders", type=Path, default=DEFAULT_REFSET / "CANONICAL_LADDER.tsv")
    parser.add_argument("--folds", type=Path, default=DEFAULT_REFSET / "DNA02_02_OUTER_FOLDS_v0.1.0.tsv")
    parser.add_argument("--sequence-objects", type=Path, default=DEFAULT_REFSET / "SEQUENCE_OBJECT.tsv")
    parser.add_argument("--fasta", type=Path, default=DEFAULT_REFSET / "DNA02_02_SEQUENCE_OBJECTS.fasta")
    parser.add_argument("--controls", type=Path, default=DEFAULT_REFSET / "DNA02_02_TAU_CONTROLS_v0.1.1.tsv")
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--repeat-embeddings", type=Path, required=True)
    parser.add_argument("--embedding-index", type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("output directory must be empty")

    formal_ids = validate_sequence_universe(args.sequence_objects, args.fasta)
    ladders = read_tsv(args.ladders)
    folds = read_tsv(args.folds)
    controls = read_tsv(args.controls)
    cluster_fold: dict[str, int] = {}
    for row in folds:
        cluster = row["outer_cluster_id"]
        fold = int(row["outer_fold"])
        if cluster in cluster_fold or fold not in range(5):
            raise ValueError("invalid frozen cluster fold assignment")
        cluster_fold[cluster] = fold

    ladder_ids: set[str] = set()
    query_meta: dict[str, tuple[str, str]] = {}
    required_objects: set[str] = set()
    for row in ladders:
        if row["ladder_id"] in ladder_ids:
            raise ValueError("duplicate ladder ID")
        ladder_ids.add(row["ladder_id"])
        if row["outer_cluster_id"] not in cluster_fold:
            raise ValueError("ladder cluster absent from frozen folds")
        if (row["hamming_x0_x1"], row["hamming_x1_x2"], row["hamming_x0_x2"]) != ("1", "1", "2"):
            raise ValueError("invalid Hamming truth")
        meta = (row["parent_locus_id"], row["outer_cluster_id"])
        if query_meta.setdefault(row["first_edit_event_id"], meta) != meta:
            raise ValueError("first-edit query identity conflict")
        required_objects.update(row[key] for key in ("x0_object_id", "x1_object_id", "x2_object_id"))

    primary = load_embeddings(args.embeddings, args.embedding_index)
    repeat = load_embeddings(args.repeat_embeddings, args.embedding_index)
    if set(primary) != formal_ids or set(primary) != required_objects or set(repeat) != formal_ids:
        raise ValueError("embedding object universe does not equal frozen RefSet universe")
    primary_shapes = {vector.shape for vector in primary.values()}
    repeat_shapes = {vector.shape for vector in repeat.values()}
    if len(primary_shapes) != 1 or primary_shapes != repeat_shapes:
        raise ValueError("primary/repeat representation dimension mismatch")
    if (len(formal_ids), len(ladders), len(query_meta), len(cluster_fold)) != (54020, 46720, 5840, 278):
        raise ValueError("DNA02-02 frozen population count mismatch")
    tau, tau_evidence = calibrate_tau(controls, ladders, primary, repeat)

    scores = [score_ladder(row, primary, tau) for row in ladders]
    formal_f, query_scores, locus_scores, cluster_scores = aggregate_hierarchy(scores)
    bootstrap_values = bootstrap_cluster_distribution(cluster_scores, draws=20000, seed=20261014)
    low, high = np.quantile(bootstrap_values, [0.025, 0.975])
    result = {
        "status": "PASS",
        "evaluation_class": "FORMAL_PROTOCOL_EVALUATION",
        "formal_scoring_protocol": True,
        "registered_paper_result": False,
        "constraint_id": "C-DNA-02-02",
        "model_id": args.model_id,
        "formal_interface": "F_ONLY",
        "score_transform_version": "DNA02-02_v0.1.5_F_ONLY",
        "F": formal_f,
        "F_95_CI": [low, high],
        "tau_delta": tau,
        "counts": {
            "objects": len(formal_ids),
            "ladders": len(scores),
            "first_edit_queries": len(query_scores),
            "outer_clusters": len(cluster_scores),
        },
        "coverage": {
            "expected_objects": 54020,
            "scored_objects": len(formal_ids),
            "expected_ladders": 46720,
            "scored_ladders": len(scores),
            "unscored_ladders": 0,
        },
        "numeric_method": {
            "operator": "SAFE_COSINE_SIMILARITY",
            "aggregation": "ladder_to_first_edit_to_locus_to_outer_cluster",
            "bootstrap": {"draws": 20000, "seed": 20261014, "unit": "outer_cluster", "interval": "percentile_95"},
            "tau_rule": "frozen_stratified_repeat_calibration_v0.1.1",
        },
        "input_sha256": {
            "ladders": sha256(args.ladders),
            "folds": sha256(args.folds),
            "sequence_objects": sha256(args.sequence_objects),
            "fasta": sha256(args.fasta),
            "controls": sha256(args.controls),
            "embeddings": sha256(args.embeddings),
            "repeat_embeddings": sha256(args.repeat_embeddings),
        },
    }
    args.output.mkdir(parents=True, exist_ok=True)
    ladder_rows = [row.__dict__ for row in scores]
    write_tsv(args.output / "INDIVIDUAL_SCORES.tsv", ladder_rows, list(ladder_rows[0]))
    write_tsv(
        args.output / "QUERY_SCORE_LEDGER.tsv",
        [
            {
                "first_edit_event_id": key,
                "parent_locus_id": query_meta[key][0],
                "outer_cluster_id": query_meta[key][1],
                "query_F": value,
            }
            for key, value in sorted(query_scores.items())
        ],
        ["first_edit_event_id", "parent_locus_id", "outer_cluster_id", "query_F"],
    )
    np.save(args.output / "F_BOOTSTRAP.npy", bootstrap_values, allow_pickle=False)
    write_tsv(
        args.output / "LOCUS_SCORE_LEDGER.tsv",
        [
            {
                "parent_locus_id": parent,
                "outer_cluster_id": next(cluster for locus, cluster in query_meta.values() if locus == parent),
                "locus_F": value,
            }
            for parent, value in sorted(locus_scores.items())
        ],
        ["parent_locus_id", "outer_cluster_id", "locus_F"],
    )
    write_tsv(
        args.output / "CLUSTER_SCORE_LEDGER.tsv",
        [
            {"outer_cluster_id": cluster, "cluster_F": value}
            for cluster, value in sorted(cluster_scores.items())
        ],
        ["outer_cluster_id", "cluster_F"],
    )
    (args.output / "TAU_EVIDENCE.json").write_text(json.dumps(tau_evidence, indent=2, sort_keys=True) + "\n")
    (args.output / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    output_files = [
        "INDIVIDUAL_SCORES.tsv",
        "QUERY_SCORE_LEDGER.tsv",
        "LOCUS_SCORE_LEDGER.tsv",
        "CLUSTER_SCORE_LEDGER.tsv",
        "F_BOOTSTRAP.npy",
        "TAU_EVIDENCE.json",
        "FORMAL_RESULT.json",
    ]
    manifest = {
        "status": "PASS",
        "files": [
            {"path": name, "bytes": (args.output / name).stat().st_size, "sha256": sha256(args.output / name)}
            for name in output_files
        ],
    }
    (args.output / "OUTPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (args.output / "RUN_COMPLETE").write_text("PASS\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
