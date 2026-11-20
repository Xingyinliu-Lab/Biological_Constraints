"""Shared DNA06 codon-local F-only implementation used by the public scorer."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np

from dna06_numeric_core import (
    bootstrap_distribution,
    fold_map,
    hierarchical,
    inner_select,
    pair_credit_ledger,
    parent_scores,
    predict,
    ridge_logistic_fit,
    site_scores,
    site_scores_from_pair_ledger,
    standardize_fit,
)

EXPECTED = {"sites": 3920, "variants": 11506, "parents": 560, "clusters": 554}


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


def load_embeddings(path: Path) -> tuple[dict[str, int], np.ndarray]:
    archive = np.load(path, allow_pickle=False)
    if set(archive.files) != {"row_ids", "embeddings"}:
        raise ValueError("codon-local NPZ must contain exactly row_ids and embeddings")
    identifiers = [str(value) for value in archive["row_ids"]]
    raw = archive["embeddings"]
    if (
        not np.issubdtype(raw.dtype, np.number)
        or np.issubdtype(raw.dtype, np.complexfloating)
        or np.issubdtype(raw.dtype, np.bool_)
    ):
        raise ValueError("embeddings must use a real numeric dtype")
    matrix = np.asarray(raw, dtype=np.float64)
    if (
        len(identifiers) != len(set(identifiers))
        or matrix.ndim != 2
        or matrix.shape[1] <= 0
        or matrix.shape[0] != len(identifiers)
        or not np.isfinite(matrix).all()
    ):
        raise ValueError("invalid embedding rows")
    return {value: index for index, value in enumerate(identifiers)}, matrix


def score(
    refset_dir: Path,
    embeddings_path: Path,
    model_id: str,
    output_dir: Path,
    evaluation_scope: str,
) -> dict:
    """Run the frozen DNA06 scientific calculation after caller-level gating."""
    refset_dir = Path(refset_dir)
    embeddings_path = Path(embeddings_path)
    output_dir = Path(output_dir)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    sites = read_tsv(refset_dir / "codon_site.tsv")
    variants = read_tsv(refset_dir / "variant_endpoint.tsv")
    members = read_tsv(refset_dir / "site_variant_membership.tsv")
    parents = read_tsv(refset_dir / "parent_cds.tsv")
    objects = read_tsv(refset_dir / "sequence_object.tsv")

    clusters = {row["parent_cds_id"]: row["dependency_cluster_id"] for row in parents}
    wt_object = {row["sequence_sha256"]: row["object_id"] for row in objects}
    parent_wt = {row["parent_cds_id"]: wt_object[row["cds_sequence_sha256"]] for row in parents}
    if (len(sites), len(variants), len(parents), len(set(clusters.values()))) != (
        EXPECTED["sites"], EXPECTED["variants"], EXPECTED["parents"], EXPECTED["clusters"]
    ):
        raise ValueError("formal population mismatch")
    if (
        len({row["site_id"] for row in sites}) != len(sites)
        or len({row["variant_id"] for row in variants}) != len(variants)
        or len({row["membership_id"] for row in members}) != len(members)
    ):
        raise ValueError("duplicate identity")

    allowed = {"SYNONYMOUS", "ORDINARY_MISSENSE"}
    site_by = {row["site_id"]: row for row in sites}
    variant_by = {row["variant_id"]: row for row in variants}
    member_by = {row["variant_id"]: row for row in members}
    if set(variant_by) != set(member_by):
        raise ValueError("incomplete ALT membership")
    for variant_id, variant in variant_by.items():
        member = member_by[variant_id]
        site = site_by.get(variant["site_id"])
        if (
            site is None
            or variant["consequence_label"] not in allowed
            or member["consequence_label"] != variant["consequence_label"]
            or member["site_id"] != variant["site_id"]
            or member["endpoint_role"] != "ALT"
            or variant["parent_cds_id"] != site["parent_cds_id"]
            or variant["changed_position_count"] != "1"
        ):
            raise ValueError("invalid matched-site bundle")

    expected_counts = {
        site_id: (int(site["synonymous_count"]), int(site["ordinary_missense_count"]))
        for site_id, site in site_by.items()
    }
    observed = {site_id: [0, 0] for site_id in site_by}
    for variant in variants:
        observed[variant["site_id"]][0 if variant["consequence_label"] == "SYNONYMOUS" else 1] += 1
    if any(tuple(observed[site_id]) != expected_counts[site_id] for site_id in site_by):
        raise ValueError("incomplete matched-site ALT set")

    index, matrix = load_embeddings(embeddings_path)
    expected_rows: set[str] = set()
    for variant in variants:
        site = site_by[variant["site_id"]]
        parent_id = variant["parent_cds_id"]
        expected_rows.add(f"{site['site_id']}|WT|{parent_wt[parent_id]}")
        expected_rows.add(
            f"{site['site_id']}|ALT|{variant['variant_id']}|{variant['object_id']}"
        )
    if set(index) != expected_rows:
        raise ValueError("missing or extra endpoint embedding rows")

    outer = fold_map(clusters.values(), 5, "DNA06_OUTER_V010")
    stage = Path(tempfile.mkdtemp(prefix=output_dir.name + ".stage.", dir=str(output_dir.parent)))
    try:
        features, labels, cluster_ids, site_ids, parent_ids, variant_ids = [], [], [], [], [], []
        for variant in variants:
            site = site_by[variant["site_id"]]
            parent_id = variant["parent_cds_id"]
            wt_key = f"{site['site_id']}|WT|{parent_wt[parent_id]}"
            alt_key = f"{site['site_id']}|ALT|{variant['variant_id']}|{variant['object_id']}"
            features.append(np.abs(matrix[index[alt_key]] - matrix[index[wt_key]]))
            labels.append(1 if variant["consequence_label"] == "ORDINARY_MISSENSE" else 0)
            cluster_ids.append(clusters[parent_id])
            site_ids.append(site["site_id"])
            parent_ids.append(parent_id)
            variant_ids.append(variant["variant_id"])

        X = np.asarray(features)
        y = np.asarray(labels)
        cluster_array = np.asarray(cluster_ids)
        site_array = np.asarray(site_ids)
        parent_array = np.asarray(parent_ids)
        predictions = np.empty(len(y))
        selected: dict[str, float] = {}
        for fold in range(5):
            train = np.array([outer[cluster] != fold for cluster in cluster_array])
            valid = ~train
            lam, _ = inner_select(
                X[train], y[train], cluster_array[train], site_array[train], parent_array[train]
            )
            mean, scale = standardize_fit(X[train])
            beta = ridge_logistic_fit((X[train] - mean) / scale, y[train], lam)
            predictions[valid] = predict((X[valid] - mean) / scale, beta)
            selected[str(fold)] = lam
        if not np.isfinite(X).all() or not np.isfinite(predictions).all():
            raise FloatingPointError("non-finite formal values")

        prediction_rows = [
            {
                "variant_id": variant_ids[index_value],
                "site_id": site_array[index_value],
                "parent_cds_id": parent_array[index_value],
                "dependency_cluster_id": cluster_array[index_value],
                "outer_fold": outer[cluster_array[index_value]],
                "label": "ORDINARY_MISSENSE" if y[index_value] else "SYNONYMOUS",
                "decision_score": format(predictions[index_value], ".17g"),
            }
            for index_value in range(len(y))
        ]
        scoring_rows = [dict(row, score=row["decision_score"]) for row in prediction_rows]
        pair_rows = pair_credit_ledger(scoring_rows)
        site_rows = site_scores_from_pair_ledger(pair_rows)
        if site_rows != site_scores(scoring_rows):
            raise AssertionError("pair-to-site aggregation mismatch")
        parent_rows = parent_scores(site_rows)
        F, cluster_scores = hierarchical(site_rows)
        by_cluster = {}
        for row in parent_rows:
            by_cluster.setdefault(row["dependency_cluster_id"], []).append(float(row["parent_F"]))
        expected_F = sum(sum(values) / len(values) for values in by_cluster.values()) / len(by_cluster)
        if abs(F - expected_F) > 1e-12:
            raise AssertionError("parent-to-cluster aggregation mismatch")
        F_bootstrap = bootstrap_distribution(cluster_scores)
        f_ci = [float(np.quantile(F_bootstrap, .025)), float(np.quantile(F_bootstrap, .975))]
        if (
            not all(np.isfinite([F, *f_ci]))
            or not -1 <= F <= 1
        ):
            raise FloatingPointError("invalid formal score")

        write_tsv(stage / "INDIVIDUAL_SCORES.tsv", prediction_rows, list(prediction_rows[0]))
        write_tsv(stage / "PAIR_CREDIT_LEDGER.tsv", pair_rows, list(pair_rows[0]))
        write_tsv(stage / "SITE_SCORE_LEDGER.tsv", site_rows, list(site_rows[0]))
        write_tsv(stage / "PARENT_SCORE_LEDGER.tsv", parent_rows, list(parent_rows[0]))
        cluster_rows = [
            {"dependency_cluster_id": key, "cluster_F": format(value, ".17g")}
            for key, value in sorted(cluster_scores.items())
        ]
        write_tsv(stage / "CLUSTER_SCORE_LEDGER.tsv", cluster_rows, list(cluster_rows[0]))

        result = {
            "constraint_id": "DNA06",
            "method_version": "0.1.1",
            "model_id": model_id,
            "evaluation_scope": evaluation_scope,
            "formal_scoring_protocol": True,
            "registered_paper_result": False,
            "score_status": "REPORTED",
            "formal_readout": "CODON_LOCAL_NUCLEOTIDE_COORDINATE",
            "representation_status": "HELD_OUT_L2_LOGISTIC_RIDGE_CHARACTERIZATION",
            "F": F,
            "F_95_CI": f_ci,
            "formal_interface": "F_ONLY",
            "score_transform_version": "DNA06_v0.1.3_F_ONLY",
            "selected_lambda_by_outer_fold": selected,
            "population": EXPECTED,
            "coverage": {
                "planned_variants": 11506,
                "scored_variants": 11506,
                "failed_variants": 0,
                "planned_sites": 3920,
                "scored_sites": 3920,
                "failed_sites": 0,
                "planned_parents": EXPECTED["parents"],
                "scored_parents": len(parent_rows),
                "scored_matched_pairs": len(pair_rows),
            },
            "embedding_sha256": sha256(embeddings_path),
        }
        result_name = "FORMAL_RESULT.json"
        (stage / result_name).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        np.save(stage / "F_BOOTSTRAP.npy", F_bootstrap, allow_pickle=False)
        files = {
            path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in stage.iterdir()
        }
        (stage / "OUTPUT_MANIFEST.json").write_text(
            json.dumps(files, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (stage / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")
        os.replace(stage, output_dir)
        return result
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
