from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scoring"))
from dna0202_numeric_core import LadderScore, aggregate_hierarchy, bootstrap_cluster_distribution
from dna0202_io import write_tsv


FROZEN_LADDER_SHA256 = "f27a6d40284997146631ed39414cae05e4872114ac347ba69f8ed9819d8386f9"
FROZEN_FASTA_SHA256 = "25fa18aed50409ba0c923ebcb0415edd2f2c93d4832263f8a6fba6c01bb3bd62"
DEFAULT_REFSET = Path(__file__).resolve().parents[2] / "refset"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    key = None
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in records:
                    raise ValueError("duplicate FASTA id")
                records[key] = ""
            elif line:
                if key is None:
                    raise ValueError("sequence before FASTA header")
                records[key] += line.upper()
    return records


def hamming(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal lengths")
    return sum(a != b for a, b in zip(left, right))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ladders", type=Path, default=DEFAULT_REFSET / "CANONICAL_LADDER.tsv")
    parser.add_argument("--fasta", type=Path, default=DEFAULT_REFSET / "DNA02_02_SEQUENCE_OBJECTS.fasta")
    parser.add_argument("--output", type=Path, required=True, help="Empty output directory")
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("output directory must be empty")
    observed_ladder_sha256 = sha256(args.ladders)
    observed_fasta_sha256 = sha256(args.fasta)
    if observed_ladder_sha256 != FROZEN_LADDER_SHA256:
        raise ValueError("ladder input does not match the frozen RefSet")
    if observed_fasta_sha256 != FROZEN_FASTA_SHA256:
        raise ValueError("FASTA input does not match the frozen RefSet")
    sequences = read_fasta(args.fasta)
    individual: list[dict[str, object]] = []
    scores: list[LadderScore] = []
    rows = 0
    with args.ladders.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            x0, x1, x2 = (sequences[row[key]] for key in ["x0_object_id", "x1_object_id", "x2_object_id"])
            observed = (hamming(x0, x1), hamming(x1, x2), hamming(x0, x2))
            if observed != (1, 1, 2):
                raise ValueError(f"invalid Level-0 truth for {row['ladder_id']}: {observed}")
            individual.append(
                {
                    "ladder_id": row["ladder_id"],
                    "first_edit_event_id": row["first_edit_event_id"],
                    "parent_locus_id": row["parent_locus_id"],
                    "outer_cluster_id": row["outer_cluster_id"],
                    "hamming_x0_x1": observed[0],
                    "hamming_x1_x2": observed[1],
                    "hamming_x0_x2": observed[2],
                    "F_credit_left": 1.0,
                    "F_credit_right": 1.0,
                    "ladder_F": 1.0,
                }
            )
            scores.append(LadderScore(
                ladder_id=row["ladder_id"], first_edit_event_id=row["first_edit_event_id"],
                parent_locus_id=row["parent_locus_id"], outer_cluster_id=row["outer_cluster_id"],
                cosine_x0_x1=1.0, cosine_x0_x2=0.0, cosine_x2_x1=1.0, cosine_x2_x0=0.0,
                delta_left=1.0, delta_right=1.0, F_credit_left=1.0, F_credit_right=1.0, ladder_F=1.0,
            ))
            rows += 1
    if rows != 46720:
        raise ValueError("formal comparator population mismatch")
    formal_f, query_scores, locus_scores, cluster_scores = aggregate_hierarchy(scores)
    query_meta = {score.first_edit_event_id: (score.parent_locus_id, score.outer_cluster_id) for score in scores}
    locus_cluster = {score.parent_locus_id: score.outer_cluster_id for score in scores}
    if formal_f != 1.0 or len(query_scores) != 5840 or len(locus_scores) != 280 or len(cluster_scores) != 278:
        raise ValueError("formal comparator hierarchy mismatch")
    bootstrap_values = bootstrap_cluster_distribution(cluster_scores, draws=20000, seed=20261014)
    low, high = np.quantile(bootstrap_values, [0.025, 0.975])
    result = {
        "status": "PASS",
        "score_status": "REPORTED",
        "constraint_id": "C-DNA-02-02",
        "model_key": "hamming_level0",
        "comparator_class": "LEVEL_0_GEOMETRIC_TRUTH",
        "formal_interface": "F_ONLY",
        "score_transform_version": "DNA02-02_v0.1.5_F_ONLY",
        "F": 1.0,
        "F_95_CI": [low, high],
        "counts": {"ladders": rows, "first_edit_queries": len(query_scores), "loci": len(locus_scores), "outer_clusters": len(cluster_scores)},
        "input_sha256": {"ladders": observed_ladder_sha256, "fasta": observed_fasta_sha256},
    }
    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output / "INDIVIDUAL_SCORES.tsv", individual, list(individual[0]))
    write_tsv(
        args.output / "QUERY_SCORE_LEDGER.tsv",
        [
            {
                "first_edit_event_id": query_id,
                "parent_locus_id": query_meta[query_id][0],
                "outer_cluster_id": query_meta[query_id][1],
                "query_F": value,
            }
            for query_id, value in sorted(query_scores.items())
        ],
        ["first_edit_event_id", "parent_locus_id", "outer_cluster_id", "query_F"],
    )
    np.save(args.output / "F_BOOTSTRAP.npy", bootstrap_values, allow_pickle=False)
    write_tsv(
        args.output / "LOCUS_SCORE_LEDGER.tsv",
        [
            {"parent_locus_id": parent, "outer_cluster_id": locus_cluster[parent], "locus_F": value}
            for parent, value in sorted(locus_scores.items())
        ],
        ["parent_locus_id", "outer_cluster_id", "locus_F"],
    )
    write_tsv(
        args.output / "CLUSTER_SCORE_LEDGER.tsv",
        [
            {"outer_cluster_id": cluster_id, "cluster_F": value}
            for cluster_id, value in sorted(cluster_scores.items())
        ],
        ["outer_cluster_id", "cluster_F"],
    )
    (args.output / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_names = ["INDIVIDUAL_SCORES.tsv", "QUERY_SCORE_LEDGER.tsv", "LOCUS_SCORE_LEDGER.tsv", "CLUSTER_SCORE_LEDGER.tsv", "F_BOOTSTRAP.npy", "FORMAL_RESULT.json"]
    manifest = {
        "status": "PASS",
        "files": [
            {"path": name, "bytes": (args.output / name).stat().st_size, "sha256": sha256(args.output / name)}
            for name in output_names
        ],
    }
    (args.output / "OUTPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.output / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "F": 1.0, "F_95_CI": [low, high]}))


if __name__ == "__main__":
    main()
