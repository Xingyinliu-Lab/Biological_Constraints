"""Run the public frozen ST04 protocol on residue representations."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import st04_numeric_core_v0_1_1 as numeric_core
from st04_io_qc import (
    constant_representation_qc,
    coverage_report,
    endpoint_set,
    load_embeddings,
    load_json,
    load_tsv,
    release_qc,
    save_fit_index,
)


def write_relation_score_ledger(rows, predictions, queries, tau, out_path):
    """Write one auditable score record for every frozen relation."""
    fields = [
        "relation_id", "query_id", "within_query_rank", "source_dependency_cluster_id",
        "structure_object_id", "anchor_residue_id", "true_target_residue_id",
        "negative_residue_id", "positive_score", "negative_score", "delta", "tie",
        "F_credit", "query_F",
    ]
    by_query_rank = {(row["query_id"], int(row["within_query_rank"])): row for row in rows}
    if len(by_query_rank) != len(rows):
        raise ValueError("ST04 relation identity closure failure")
    relation_credit = {}
    with gzip.open(out_path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for query_id in sorted(predictions):
            prediction = predictions[query_id]
            negatives = queries[query_id]["negatives"]
            if len(negatives) != len(prediction["negative_scores"]) or len(negatives) != len(prediction["credits"]):
                raise ValueError("ST04 OOF negative-score closure failure")
            credits = []
            for rank, (negative_id, negative_score, credit) in enumerate(
                zip(negatives, prediction["negative_scores"], prediction["credits"]), 1
            ):
                row = by_query_rank[(query_id, rank)]
                if row["negative_residue_id"] != negative_id:
                    raise ValueError("ST04 OOF relation endpoint closure failure")
                delta = float(prediction["positive_score"] - negative_score)
                expected_credit = 1.0 if delta > tau else -1.0 if delta < -tau else 0.0
                if float(credit) != expected_credit:
                    raise ValueError("ST04 relation credit consistency failure")
                tie = abs(delta) <= tau
                writer.writerow({
                    "relation_id": row["relation_id"], "query_id": query_id,
                    "within_query_rank": rank, "source_dependency_cluster_id": row["source_dependency_cluster_id"],
                    "structure_object_id": row["structure_object_id"], "anchor_residue_id": row["anchor_residue_id"],
                    "true_target_residue_id": row["true_target_residue_id"], "negative_residue_id": negative_id,
                    "positive_score": prediction["positive_score"], "negative_score": negative_score,
                    "delta": delta, "tie": str(tie).lower(), "F_credit": credit,
                    "query_F": prediction["query_credit"],
                })
                credits.append(float(credit))
                relation_credit[row["relation_id"]] = float(credit)
            if not np.isclose(float(np.mean(credits)), prediction["query_credit"], rtol=0.0, atol=1e-12):
                raise ValueError("ST04 query-to-relation aggregation consistency failure")
    return {
        "status": "PASS", "relations": len(relation_credit), "queries": len(predictions),
        "relation_credit_sha256": hashlib.sha256(
            "\n".join(f"{key}\t{relation_credit[key]:.17g}" for key in sorted(relation_credit)).encode()
        ).hexdigest(),
    }


def validate_relation_ledger_aggregation(ledger_path, expected_hierarchy, expected_r):
    """Reconstruct the frozen macro hierarchy from the emitted relation ledger."""
    query_rows = {}
    with gzip.open(ledger_path, "rt", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query_rows.setdefault(row["query_id"], []).append(row)
    query_scores = {}
    query_meta = {}
    for query_id, rows in query_rows.items():
        values = np.asarray([float(row["F_credit"]) for row in rows], dtype=np.float64)
        recorded = {float(row["query_F"]) for row in rows}
        if len(recorded) != 1 or not np.isclose(values.mean(), recorded.pop(), rtol=0.0, atol=1e-12):
            raise ValueError("ST04 ledger relation-to-query aggregation failure")
        first = rows[0]
        query_scores[query_id] = float(values.mean())
        query_meta[query_id] = (
            first["source_dependency_cluster_id"], first["structure_object_id"], first["anchor_residue_id"],
        )
    anchor_queries = defaultdict(list)
    for query_id, meta in query_meta.items():
        anchor_queries[meta].append(query_scores[query_id])
    anchors = {"|".join(key): float(np.mean(values)) for key, values in anchor_queries.items()}
    structure_anchors = defaultdict(list)
    for key, value in anchors.items():
        cluster, structure, _ = key.split("|", 2)
        structure_anchors[(cluster, structure)].append(value)
    structures = {"|".join(key): float(np.mean(values)) for key, values in structure_anchors.items()}
    cluster_structures = defaultdict(list)
    for key, value in structures.items():
        cluster, _ = key.split("|", 1)
        cluster_structures[cluster].append(value)
    cluster_values = {cluster: float(np.mean(values)) for cluster, values in cluster_structures.items()}
    reconstructed_r = float(np.mean(list(cluster_values.values())))
    for key in ("query", "anchor", "structure", "cluster"):
        observed = {"query": query_scores, "anchor": anchors, "structure": structures, "cluster": cluster_values}[key]
        expected = expected_hierarchy[key]
        if set(observed) != set(expected) or not all(np.isclose(observed[item], expected[item], rtol=0.0, atol=1e-12) for item in observed):
            raise ValueError(f"ST04 ledger {key} hierarchy consistency failure")
    if not np.isclose(reconstructed_r, expected_r, rtol=0.0, atol=1e-12):
        raise ValueError("ST04 ledger cluster-to-R aggregation failure")
    return {
        "status": "PASS", "queries": len(query_scores), "anchors": len(anchors),
        "structures": len(structures), "clusters": len(cluster_values), "F": reconstructed_r,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", type=Path, required=True)
    parser.add_argument("--refset-release", type=Path, required=True)
    parser.add_argument("--numeric-method", type=Path, required=True)
    parser.add_argument("--split-release", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--scope-record", type=Path,
                        help="required only for a documented model-input limitation relation table")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    method = load_json(args.numeric_method)
    refset = load_json(args.refset_release)
    split = load_json(args.split_release)
    if method.get("status") != "FROZEN" or refset.get("status") != "FROZEN" or split.get("status") != "FROZEN":
        raise ValueError("ST04 frozen method, RefSet, or split status failure")
    rows = load_tsv(args.relations)
    queries = numeric_core.build_queries(rows)
    clusters = {query["cluster"] for query in queries.values()}
    endpoints = endpoint_set(rows)
    relation_sha256 = hashlib.sha256(args.relations.read_bytes()).hexdigest()
    if relation_sha256 == refset["scientific_tables"]["relations.tsv"]:
        scope_identity = "FULL_PUBLIC_REFSET"
        expected_relations = refset["formal_population"]["formal_relations"]
        expected_queries = refset["formal_population"]["formal_directed_queries"]
        expected_clusters = refset["formal_population"]["formal_dependency_clusters"]
        expected_endpoints = 48282
    else:
        if args.scope_record is None:
            raise ValueError("ST04 non-FULL relation input requires a documented --scope-record")
        scope = load_json(args.scope_record)
        table = scope.get("relation_table", {})
        if (scope.get("constraint_id") != "ST04" or
                scope.get("scope_kind") != "MODEL_INPUT_LIMITATION" or
                table.get("sha256") != relation_sha256):
            raise ValueError("ST04 model-input limitation identity mismatch")
        scope_identity = "MODEL_INPUT_LIMITATION"
        expected_relations = table.get("relations")
        expected_queries = table.get("directed_queries")
        expected_clusters = table.get("dependency_clusters")
        expected_endpoints = table.get("residue_endpoints")
    expected = (expected_relations, expected_queries, expected_clusters, expected_endpoints)
    observed = (len(rows), len(queries), len(clusters), len(endpoints))
    if observed != expected:
        raise ValueError(f"ST04 frozen population mismatch: {observed} != {expected}")

    embeddings, dimension = load_embeddings(args.embeddings)
    if set(embeddings) != set(endpoints):
        raise ValueError("embedding residue endpoint closure failure")
    if set(split["outer"]) != clusters:
        raise ValueError("outer split cluster closure failure")
    for outer in range(method["cross_fitting"]["outer_folds"]):
        training_clusters = {cluster for cluster in clusters if split["outer"][cluster] != outer}
        if set(split["inner"][str(outer)]) != training_clusters:
            raise ValueError("inner split cluster closure failure")

    qc = release_qc(rows, queries, embeddings)
    if qc["status"] != "PASS":
        raise ValueError("ST04 representation or RefSet QC failure")
    numeric_method = {
        "outer_folds": method["cross_fitting"]["outer_folds"],
        "inner_folds": method["cross_fitting"]["inner_folds"],
        "lambda_grid": method["training"]["lambda_grid"],
        "solver": method["training"]["solver"],
        "tie_tolerance": method["held_out_credit"]["tie_tolerance"],
    }
    formal_f, predictions, cluster_scores, selection, fits, hierarchy = numeric_core.nested_crossfit(
        queries, embeddings, split, numeric_method
    )
    ci = numeric_core.bootstrap_cluster_means(
        cluster_scores, method["uncertainty"]["draws"], method["uncertainty"]["seed"]
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fit_index = save_fit_index(fits, args.out_dir)
    (args.out_dir / "RELEASE_QC.json").write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
    (args.out_dir / "SCOPE_AND_STRATIFIED_COVERAGE.json").write_text(
        json.dumps(coverage_report(rows, queries), indent=2, sort_keys=True) + "\n"
    )
    (args.out_dir / "LAMBDA_SELECTION.json").write_text(json.dumps(selection, indent=2, sort_keys=True) + "\n")
    (args.out_dir / "HIERARCHY_SCORES.json").write_text(json.dumps(hierarchy, indent=2, sort_keys=True) + "\n")
    with (args.out_dir / "OOF_QUERY_SCORES.jsonl").open("w") as handle:
        for query_id in sorted(predictions):
            handle.write(json.dumps({"query_id": query_id, **predictions[query_id]}, sort_keys=True) + "\n")
    relation_ledger = write_relation_score_ledger(
        rows, predictions, queries, numeric_method["tie_tolerance"], args.out_dir / "RELATION_SCORE_LEDGER.tsv.gz"
    )
    relation_aggregation = validate_relation_ledger_aggregation(
        args.out_dir / "RELATION_SCORE_LEDGER.tsv.gz", hierarchy, formal_f
    )
    (args.out_dir / "RELATION_LEDGER_VALIDATION.json").write_text(
        json.dumps({"relation": relation_ledger, "aggregation": relation_aggregation}, indent=2, sort_keys=True) + "\n"
    )
    input_sha256 = {
        "embeddings": hashlib.sha256(args.embeddings.read_bytes()).hexdigest(),
        "relations": hashlib.sha256(args.relations.read_bytes()).hexdigest(),
        "refset_release": hashlib.sha256(args.refset_release.read_bytes()).hexdigest(),
        "numeric_method": hashlib.sha256(args.numeric_method.read_bytes()).hexdigest(),
        "split_release": hashlib.sha256(args.split_release.read_bytes()).hexdigest(),
    }
    result = {
        "status": "PASS",
        "evaluation_class": "FORMAL_PROTOCOL_EVALUATION",
        "formal_scoring_protocol": True,
        "registered_paper_result": False,
        "scope_identity": scope_identity,
        "constraint_id": "C-ST-04",
        "model_id": args.model_id,
        "formal_interface": "F_ONLY",
        "score_transform_version": "ST04_v0.1.4_F_ONLY",
        "F": formal_f,
        "F_95_CI": ci,
        "counts": {
            "relations": len(rows),
            "directed_queries": len(queries),
            "dependency_clusters": len(cluster_scores),
            "required_residue_endpoints": len(endpoints),
            "embedding_dimension": dimension,
        },
        "fit_parameters_sha256": fit_index["fit_parameters_sha256"],
        "input_sha256": input_sha256,
        "relation_score_ledger": {"relation": relation_ledger, "aggregation": relation_aggregation},
    }
    (args.out_dir / "EVALUATION_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
