#!/usr/bin/env python3
"""Generate the DNA03 Level-0 reference through the public aggregation core."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import numpy as np
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

CONSTRAINT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CONSTRAINT_ROOT / "scoring"))

from dna03_numeric_core import aggregate_query_rows, cluster_bootstrap_distribution  # noqa: E402

EXPECTED = {"queries": 1120, "pairs": 35840, "dependency_clusters": 278, "transitions": 4}
DEFAULT_REFSET = CONSTRAINT_ROOT / "refset"
BOOTSTRAP_DRAWS = 20000
BOOTSTRAP_SEED = 20261014


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    record_id: str | None = None
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if record_id is not None:
                records[record_id] = "".join(sequence).upper()
            record_id = line[1:].split()[0]
            sequence = []
        elif line.strip():
            sequence.append(line.strip())
    if record_id is not None:
        records[record_id] = "".join(sequence).upper()
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def containment_distance(first: str, second: str) -> int:
    short, long = (first, second) if len(first) <= len(second) else (second, first)
    return 0 if short in long or reverse_complement(short) in long else 1


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--refset", type=Path, default=DEFAULT_REFSET)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_REFSET / "DNA03_FORWARD_VIEWS.fasta")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError(f"output directory must be empty: {args.output}")

    query_rows_source = read_tsv(args.refset / "PRIMARY_QUERY.tsv")
    negatives = read_tsv(args.refset / "PRIMARY_NEGATIVE.tsv")
    if len({row["query_id"] for row in query_rows_source}) != len(query_rows_source):
        raise ValueError("duplicate query_id")
    if len({row["pair_id"] for row in negatives}) != len(negatives):
        raise ValueError("duplicate pair_id")
    negatives_by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in negatives:
        negatives_by_query[row["query_id"]].append(row)
    sequences = read_fasta(args.fasta)

    pair_rows: list[dict] = []
    query_rows: list[dict] = []
    for source in query_rows_source:
        candidates = sorted(negatives_by_query[source["query_id"]], key=lambda row: int(row["negative_rank"]))
        if len(candidates) != 32:
            raise ValueError("negative candidate closure failure")
        core = sequences[source["query_view_id"]]
        mate = sequences[source["true_mate_view_id"]]
        positive_d0 = containment_distance(core, mate)
        credits: list[float] = []
        for candidate in candidates:
            negative_d0 = containment_distance(core, sequences[candidate["negative_view_id"]])
            F_credit = 1.0 if positive_d0 < negative_d0 else 0.0 if positive_d0 == negative_d0 else -1.0
            credits.append(F_credit)
            pair_rows.append(
                {
                    "pair_id": candidate["pair_id"],
                    "query_id": source["query_id"],
                    "core_id": source["core_id"],
                    "transition_id": source["transition_id"],
                    "dependency_cluster_id": source["dependency_cluster_id"],
                    "positive_d0": positive_d0,
                    "negative_d0": negative_d0,
                    "F_credit": F_credit,
                }
            )
        query_rows.append(
            {
                "query_id": source["query_id"],
                "core_id": source["core_id"],
                "transition_id": source["transition_id"],
                "dependency_cluster_id": source["dependency_cluster_id"],
                "M": len(credits),
                "query_F": sum(credits) / len(credits),
            }
        )

    F, transition_scores, cluster_scores = aggregate_query_rows(query_rows)
    F_bootstrap = cluster_bootstrap_distribution(query_rows, BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)
    ci = [float(np.quantile(F_bootstrap, .025)), float(np.quantile(F_bootstrap, .975))]
    cluster_counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in query_rows:
        cluster_counts[(row["dependency_cluster_id"], row["transition_id"])] += 1
    cluster_rows = [
        {
            "dependency_cluster_id": cluster,
            "transition_id": transition,
            "core_count": cluster_counts[(cluster, transition)],
            "cluster_transition_F": value,
        }
        for (cluster, transition), value in sorted(cluster_scores.items())
    ]
    counts = {
        "queries": len(query_rows),
        "pairs": len(pair_rows),
        "dependency_clusters": len({row["dependency_cluster_id"] for row in query_rows}),
        "transitions": len(transition_scores),
    }
    truth = {
        "positive_containment_all": all(row["positive_d0"] == 0 for row in pair_rows),
        "negative_noncontainment_all": all(row["negative_d0"] == 1 for row in pair_rows),
        "all_pair_credits_one": all(row["F_credit"] == 1.0 for row in pair_rows),
    }
    if counts != EXPECTED or not all(truth.values()):
        raise ValueError("formal population/truth validation failed")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{args.output.name}.staging.", dir=args.output.parent))
    try:
        write_tsv(stage / "INDIVIDUAL_SCORES.tsv", pair_rows, ["pair_id", "query_id", "core_id", "transition_id", "dependency_cluster_id", "positive_d0", "negative_d0", "F_credit"])
        write_tsv(stage / "QUERY_SCORE_LEDGER.tsv", query_rows, ["query_id", "core_id", "transition_id", "dependency_cluster_id", "M", "query_F"])
        write_tsv(stage / "CLUSTER_TRANSITION_SCORE_LEDGER.tsv", cluster_rows, ["dependency_cluster_id", "transition_id", "core_count", "cluster_transition_F"])
        result = {
            "status": "PASS", "evaluation_class": "FORMAL_COMPARATOR_REFERENCE", "formal_scoring_protocol": True,
            "constraint_id": "C-DNA-03", "formal_readout": "CORE_LOCAL_NUCLEOTIDE_COORDINATE",
            "comparator_id": "EXACT_CORE_CONTAINMENT_LEVEL0", "comparator_class": "LEVEL_0_GEOMETRIC_TRUTH",
            "formal_interface": "F_ONLY", "score_transform_version": "DNA03_v0.2.6_F_ONLY",
            "F": F, "F_95_CI": ci, "transition_F": transition_scores,
            "counts": counts, "coverage": {"expected_pairs": EXPECTED["pairs"], "scored_pairs": len(pair_rows), "unscored_pairs": 0},
            "bootstrap": {"unit": "DEPENDENCY_CLUSTER", "draws": BOOTSTRAP_DRAWS, "seed": BOOTSTRAP_SEED, "ci": "PERCENTILE_95"},
            "input_sha256": {"PRIMARY_QUERY.tsv": sha256(args.refset / "PRIMARY_QUERY.tsv"), "PRIMARY_NEGATIVE.tsv": sha256(args.refset / "PRIMARY_NEGATIVE.tsv"), "DNA03_FORWARD_VIEWS.fasta": sha256(args.fasta)},
            "truth_checks": truth,
        }
        (stage / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        np.save(stage / "F_BOOTSTRAP.npy", F_bootstrap, allow_pickle=False)
        names = ["INDIVIDUAL_SCORES.tsv", "QUERY_SCORE_LEDGER.tsv", "CLUSTER_TRANSITION_SCORE_LEDGER.tsv", "F_BOOTSTRAP.npy", "FORMAL_RESULT.json"]
        manifest = {"status": "PASS", "files": [
            {"path": name, "bytes": (stage / name).stat().st_size, "sha256": sha256(stage / name)}
            for name in names
        ]}
        (stage / "OUTPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (stage / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")
        os.replace(stage, args.output)
        print(json.dumps(result, indent=2, sort_keys=True))
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
