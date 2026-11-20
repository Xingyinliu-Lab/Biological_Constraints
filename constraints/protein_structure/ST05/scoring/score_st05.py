"""Run the public frozen ST05 invariant scoring protocol."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import st05_f_only_core as numeric_core
from st05_scoring_core import evaluate_law, string_array


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def real_numeric_array(values: np.ndarray, label: str) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.kind not in {"i", "u", "f"} or array.dtype.kind == "b":
        raise ValueError(f"{label} must have a real integer or floating dtype")
    return np.asarray(array, dtype=np.float64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relations", type=Path, required=True)
    parser.add_argument("--refset-release", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output-dtype", choices=("float32", "float64"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    refset = json.loads(args.refset_release.read_text())
    if refset.get("status") != "FROZEN":
        raise ValueError("ST05 RefSet release is not frozen")
    with args.relations.open(newline="") as handle:
        relations = list(csv.DictReader(handle, delimiter="\t"))
    counts = refset["counts"]
    if len(relations) != counts["relations"]:
        raise ValueError("ST05 relation count mismatch")
    relation_ids = [row["relation_id"] for row in relations]
    if len(relation_ids) != len(set(relation_ids)):
        raise ValueError("duplicate relation ID")
    base_ids_expected = {row["base_object_id"] for row in relations}
    transform_ids = {row["transform_id"] for row in relations}
    if len(base_ids_expected) != counts["base_objects"]:
        raise ValueError("base-object closure failure")
    bundles: dict[str, set[str]] = defaultdict(set)
    for row in relations:
        bundles[row["base_object_id"]].add(row["transform_id"])
    if any(bundle != transform_ids for bundle in bundles.values()):
        raise ValueError("incomplete transform bundle")

    with np.load(args.embeddings, allow_pickle=False) as data:
        required = {"base_object_ids", "relation_ids", "base_embeddings", "variant_embeddings"}
        if set(data.files) != required:
            raise ValueError(f"embedding NPZ must contain exactly {sorted(required)}")
        base_ids = string_array(data["base_object_ids"])
        observed_relation_ids = string_array(data["relation_ids"])
        raw_base_embeddings = data["base_embeddings"]
        raw_variant_embeddings = data["variant_embeddings"]
        declared_dtype = np.dtype(args.output_dtype)
        if raw_base_embeddings.dtype != declared_dtype or raw_variant_embeddings.dtype != declared_dtype:
            raise ValueError(
                "embedding dtype must match --output-dtype before scorer conversion: "
                f"base={raw_base_embeddings.dtype}, variant={raw_variant_embeddings.dtype}, "
                f"declared={declared_dtype}"
            )
        base_embeddings = real_numeric_array(raw_base_embeddings, "base_embeddings")
        variant_embeddings = real_numeric_array(raw_variant_embeddings, "variant_embeddings")
    if len(base_ids) != len(set(base_ids)) or set(base_ids) != base_ids_expected:
        raise ValueError("embedding base-object closure failure")
    if len(observed_relation_ids) != len(set(observed_relation_ids)) or set(observed_relation_ids) != set(relation_ids):
        raise ValueError("embedding relation closure failure")
    if (
        base_embeddings.ndim != 2
        or variant_embeddings.ndim != 2
        or base_embeddings.shape != (len(base_ids), base_embeddings.shape[1])
        or variant_embeddings.shape != (len(observed_relation_ids), base_embeddings.shape[1])
        or base_embeddings.shape[1] == 0
    ):
        raise ValueError("embedding array shape failure")
    if not np.isfinite(base_embeddings).all() or not np.isfinite(variant_embeddings).all():
        raise ValueError("embeddings contain NaN or Inf")

    dimension = base_embeddings.shape[1]
    epsilon = np.finfo(np.float32 if args.output_dtype == "float32" else np.float64).eps
    tau_abs = float(256.0 * np.sqrt(dimension) * epsilon)
    tau_rel = tau_abs
    base_map = dict(zip(base_ids, base_embeddings, strict=True))
    variant_map = dict(zip(observed_relation_ids, variant_embeddings, strict=True))
    identity = np.eye(dimension, dtype=np.float64)
    zero = np.zeros(dimension, dtype=np.float64)

    scored_rows = []
    residuals = []
    ratios = []
    strata: dict[str, list[float]] = defaultdict(list)
    for relation in relations:
        evaluated = evaluate_law(
            base_map[relation["base_object_id"]],
            variant_map[relation["relation_id"]],
            identity,
            zero,
            None,
            tau_abs,
            tau_rel,
        )
        f_credit = evaluated["credit"]
        scored_rows.append(
            {
                **relation,
                "residual": evaluated["residual"],
                "scale": evaluated["scale"],
                "tolerance": evaluated["tau"],
                "normalized_residual": evaluated["ratio"],
                "margin": evaluated["tau"] - evaluated["residual"],
                "F_credit": f_credit,
                "status": "PASS" if f_credit else "LAW_VIOLATION",
                "threshold_rule": "residual <= tolerance",
            }
        )
        residuals.append(evaluated["residual"])
        ratios.append(evaluated["ratio"])
        strata[relation["sampling_stratum"]].append(f_credit)

    result = numeric_core.aggregate(scored_rows)
    ci = numeric_core.bootstrap_interval(result["cluster_scores"], draws=20000, seed=15005)
    identity_rows = [row for row in scored_rows if row["is_identity_transform"].lower() == "true"]
    if len(identity_rows) != len(base_ids_expected) or any(row["F_credit"] != 1.0 for row in identity_rows):
        raise ValueError("identity-transform QC failure")
    output = {
        "status": "PASS",
        "evaluation_class": "FORMAL_PROTOCOL_EVALUATION",
        "formal_scoring_protocol": True,
        "registered_paper_result": False,
        "constraint_id": "C-ST-05",
        "model_id": args.model_id,
        "branch": "INVARIANT",
        "formal_interface": "F_ONLY",
        "score_transform_version": "ST05_v0.1.2_F_ONLY_IDENTITY",
        "F": result["F"],
        "F_95_CI": list(ci),
        "counts": counts,
        "numerical_tolerance": {
            "output_dtype": args.output_dtype,
            "dimension": dimension,
            "formula": "256*sqrt(D)*machine_epsilon(output_dtype)",
            "tau_abs": tau_abs,
            "tau_rel": tau_rel,
        },
        "raw_residual_summary": {
            "min": float(np.min(residuals)),
            "median": float(np.median(residuals)),
            "max": float(np.max(residuals)),
        },
        "residual_tolerance_ratio_summary": {
            "min": float(np.min(ratios)),
            "median": float(np.median(ratios)),
            "max": float(np.max(ratios)),
        },
        "stratum_F": {key: float(np.mean(values)) for key, values in sorted(strata.items())},
        "hierarchy_F": {
            "base_F": result["base_scores"],
            "parent_F": result["parent_scores"],
            "dependency_cluster_F": result["cluster_scores"],
        },
        "input_identities": {
            "embeddings_sha256": sha256(args.embeddings),
            "relations_sha256": sha256(args.relations),
            "refset_release_sha256": sha256(args.refset_release),
            "scorer_sha256": sha256(Path(__file__)),
            "numeric_core_sha256": sha256(Path(numeric_core.__file__)),
            "law_core_sha256": sha256(Path(evaluate_law.__code__.co_filename)),
        },
        "identity_qc": "PASS",
        "bundle_completeness_qc": "PASS",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    ledger = args.out.with_name("ST05_TRANSFORM_RESULT_LEDGER.tsv")
    fields = list(relations[0]) + ["residual", "scale", "tolerance", "normalized_residual", "margin", "F_credit", "status", "threshold_rule"]
    with ledger.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(scored_rows)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, sort_keys=True))


if __name__ == "__main__":
    main()
