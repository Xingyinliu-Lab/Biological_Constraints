"""Run the public frozen DNA03 core-local scoring protocol."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

from dna03_numeric_core import aggregate_query_rows, cluster_bootstrap_distribution, cosine_distance, tie_F_credit
from dna03_tau_core import MAX_AUTHORIZED_TAU, compute_tau, load_embeddings

EXPECTED = {"queries": 1120, "pairs": 35840, "clusters": 278, "transitions": 4}
BOOTSTRAP_DRAWS = 20000
BOOTSTRAP_SEED = 20261014
CONSTRAINT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFSET = CONSTRAINT_ROOT / "refset"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Score core-local representations against DNA03.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--queries", default=DEFAULT_REFSET / "PRIMARY_QUERY.tsv", type=Path)
    parser.add_argument("--negatives", default=DEFAULT_REFSET / "PRIMARY_NEGATIVE.tsv", type=Path)
    parser.add_argument("--embeddings", required=True, type=Path)
    parser.add_argument("--repeat-embeddings", required=True, type=Path)
    parser.add_argument("--controls", default=DEFAULT_REFSET / "DNA03_TAU_CONTROL_IDS_v0.1.1.tsv", type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("output directory must be empty")

    tau_evidence = compute_tau(args.controls, args.embeddings, args.repeat_embeddings)
    tau = float(tau_evidence["selected_tau"])
    if tau_evidence["maximum_abs_repeat_contrast_delta"] > MAX_AUTHORIZED_TAU:
        raise ValueError("unstable repeat interface")

    queries = read_tsv(args.queries)
    negatives = read_tsv(args.negatives)
    if len({row["pair_id"] for row in negatives}) != len(negatives):
        raise ValueError("duplicate pair_id")
    negatives_by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in negatives:
        negatives_by_query[row["query_id"]].append(row)
    observed = (
        len(queries),
        len(negatives),
        len({row["dependency_cluster_id"] for row in queries}),
        len({row["transition_id"] for row in queries}),
    )
    if observed != (1120, 35840, 278, 4):
        raise ValueError("formal population mismatch")
    if len({row["query_id"] for row in queries}) != len(queries):
        raise ValueError("duplicate query_id")
    if any(len(negatives_by_query[row["query_id"]]) != 32 for row in queries):
        raise ValueError("negative candidate closure failure")

    embedding_index, arrays = load_embeddings(args.embeddings)
    repeat_index, _ = load_embeddings(args.repeat_embeddings)
    if set(embedding_index) != set(repeat_index):
        raise ValueError("primary/repeat view universes differ")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{args.output.name}.staging.", dir=args.output.parent))
    try:
        array = arrays
        pair_rows: list[dict] = []
        query_rows: list[dict] = []
        required_views: set[str] = set()
        for query_row in queries:
            candidates = negatives_by_query[query_row["query_id"]]
            required = [
                query_row["query_view_id"],
                query_row["true_mate_view_id"],
                *[candidate["negative_view_id"] for candidate in candidates],
            ]
            required_views.update(required)
            if any(view_id not in embedding_index for view_id in required):
                raise ValueError("embedding coverage mismatch")
            query_vector = array[embedding_index[query_row["query_view_id"]]]
            positive_distance = cosine_distance(
                query_vector, array[embedding_index[query_row["true_mate_view_id"]]]
            )
            credits = []
            for candidate in candidates:
                negative_distance = cosine_distance(
                    query_vector, array[embedding_index[candidate["negative_view_id"]]]
                )
                margin = negative_distance - positive_distance
                F_credit = tie_F_credit(positive_distance, negative_distance, tau)
                credits.append(F_credit)
                pair_rows.append(
                    {
                        "pair_id": candidate["pair_id"],
                        "query_id": query_row["query_id"],
                        "core_id": query_row["core_id"],
                        "transition_id": query_row["transition_id"],
                        "dependency_cluster_id": query_row["dependency_cluster_id"],
                        "query_view_id": query_row["query_view_id"],
                        "true_mate_view_id": query_row["true_mate_view_id"],
                        "negative_view_id": candidate["negative_view_id"],
                        "d_positive": positive_distance,
                        "d_negative": negative_distance,
                        "margin": margin,
                        "tau": tau,
                        "tie": abs(margin) <= tau,
                        "F_credit": F_credit,
                    }
                )
            query_rows.append(
                {
                    "query_id": query_row["query_id"],
                    "core_id": query_row["core_id"],
                    "transition_id": query_row["transition_id"],
                    "dependency_cluster_id": query_row["dependency_cluster_id"],
                    "M": len(credits),
                    "query_F": float(np.mean(credits)),
                }
            )
        if set(embedding_index) != required_views:
            raise ValueError("embedding view universe contains missing or extra views")

        pair_by_query: dict[str, list[float]] = defaultdict(list)
        for row in pair_rows:
            pair_by_query[row["query_id"]].append(float(row["F_credit"]))
        if any(abs(float(row["query_F"]) - float(np.mean(pair_by_query[row["query_id"]]))) > 1e-12 for row in query_rows):
            raise AssertionError("pair-to-query aggregation mismatch")
        F, transition_scores, cluster_scores = aggregate_query_rows(query_rows)
        F_bootstrap = cluster_bootstrap_distribution(query_rows, BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)
        ci = [float(np.quantile(F_bootstrap, .025)), float(np.quantile(F_bootstrap, .975))]
        cluster_rows = [
            {
                "dependency_cluster_id": cluster,
                "transition_id": transition,
                "cluster_transition_F": value,
            }
            for (cluster, transition), value in sorted(cluster_scores.items())
        ]
        transition_rows = [
            {"transition_id": transition, "cluster_count": sum(1 for row in cluster_rows if row["transition_id"] == transition), "transition_F": value}
            for transition, value in sorted(transition_scores.items())
        ]
        if abs(F - float(np.mean([row["transition_F"] for row in transition_rows]))) > 1e-12:
            raise AssertionError("cluster-to-transition aggregation mismatch")
        write_tsv(
            stage / "QUERY_SCORE_LEDGER.tsv",
            query_rows,
            ["query_id", "core_id", "transition_id", "dependency_cluster_id", "M", "query_F"],
        )
        write_tsv(
            stage / "INDIVIDUAL_SCORES.tsv",
            pair_rows,
            [
                "pair_id", "query_id", "core_id", "transition_id", "dependency_cluster_id",
                "query_view_id", "true_mate_view_id", "negative_view_id",
                "d_positive", "d_negative", "margin", "tau", "tie", "F_credit",
            ],
        )
        write_tsv(
            stage / "CLUSTER_TRANSITION_SCORE_LEDGER.tsv",
            cluster_rows,
            ["dependency_cluster_id", "transition_id", "cluster_transition_F"],
        )
        write_tsv(stage / "TRANSITION_SCORE_LEDGER.tsv", transition_rows, ["transition_id", "cluster_count", "transition_F"])
        result = {
            "status": "PASS",
            "evaluation_class": "FORMAL_PROTOCOL_EVALUATION",
            "formal_scoring_protocol": True,
            "registered_paper_result": False,
            "constraint_id": "C-DNA-03",
            "model_id": args.model_id,
            "formal_interface": "F_ONLY",
            "score_transform_version": "DNA03_v0.2.6_F_ONLY",
            "formal_readout": "CORE_LOCAL_NUCLEOTIDE_COORDINATE",
            "F": F,
            "F_95_CI": ci,
            "transition_F": transition_scores,
            "tau": tau,
            "counts": EXPECTED,
            "coverage": {
                "expected_views": len(required_views),
                "scored_views": len(required_views),
                "expected_queries": EXPECTED["queries"],
                "scored_queries": len(query_rows),
                "expected_pairs": EXPECTED["pairs"],
                "scored_pairs": len(pair_rows),
                "unscored_pairs": 0,
            },
            "numeric_method": {
                "distance": "1-cosine_similarity",
                "aggregation": "pair_to_query_to_dependency_cluster_transition_to_transition_macro",
                "bootstrap": {
                    "unit": "dependency_cluster",
                    "draws": BOOTSTRAP_DRAWS,
                    "seed": BOOTSTRAP_SEED,
                    "interval": "percentile_95",
                },
                "tau_rule": "frozen_stratified_repeat_calibration_v0.2.5",
            },
            "repeat_validation": {
                "controls": 128,
                "maximum_abs_repeat_contrast_delta": tau_evidence[
                    "maximum_abs_repeat_contrast_delta"
                ],
                "status": "PASS",
            },
            "input_sha256": {
                "PRIMARY_QUERY.tsv": sha256(args.queries),
                "PRIMARY_NEGATIVE.tsv": sha256(args.negatives),
                "embeddings": sha256(args.embeddings),
                "repeat_embeddings": sha256(args.repeat_embeddings),
                "controls": sha256(args.controls),
            },
        }
        (stage / "FORMAL_RESULT.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (stage / "TAU_EVIDENCE.json").write_text(
            json.dumps(tau_evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        np.save(stage / "F_BOOTSTRAP.npy", F_bootstrap, allow_pickle=False)
        output_names = [
            "INDIVIDUAL_SCORES.tsv",
            "QUERY_SCORE_LEDGER.tsv",
            "CLUSTER_TRANSITION_SCORE_LEDGER.tsv",
            "TRANSITION_SCORE_LEDGER.tsv",
            "TAU_EVIDENCE.json",
            "F_BOOTSTRAP.npy",
            "FORMAL_RESULT.json",
        ]
        manifest = {
            "status": "PASS",
            "files": [
                {"path": name, "bytes": (stage / name).stat().st_size, "sha256": sha256(stage / name)}
                for name in output_names
            ],
        }
        (stage / "OUTPUT_MANIFEST.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (stage / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")
        os.replace(stage, args.output)
        print(json.dumps(result, sort_keys=True))
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
