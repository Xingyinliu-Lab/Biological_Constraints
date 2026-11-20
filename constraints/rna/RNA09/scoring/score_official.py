#!/usr/bin/env python3
"""Score RNA09 from a validated RefSet and model representations."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

BOOTSTRAP_SEED = 20261004
BOOTSTRAP_DRAWS = 20000


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows):
    if not rows:
        raise ValueError("cannot write an empty ledger")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path):
    records, identifier, sequence = {}, None, []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if identifier in records:
                        raise ValueError(f"duplicate FASTA identifier: {identifier}")
                    records[identifier] = "".join(sequence).upper()
                identifier, sequence = line[1:].split()[0], []
                if not identifier:
                    raise ValueError("empty FASTA identifier")
            elif identifier is None:
                raise ValueError("FASTA sequence before first identifier")
            else:
                sequence.append(line)
    if identifier is not None:
        if identifier in records:
            raise ValueError(f"duplicate FASTA identifier: {identifier}")
        records[identifier] = "".join(sequence).upper()
    if not records:
        raise ValueError("empty endpoint FASTA")
    return records


def cosine(a, b):
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("non-finite cosine vector")
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    if not np.isfinite(denominator):
        raise ValueError("non-finite cosine denominator")
    if denominator == 0.0:
        return 0.0
    value = float(np.dot(a, b) / denominator)
    if not np.isfinite(value):
        raise ValueError("non-finite cosine similarity")
    return value


def assign_credits(rows):
    for row in rows:
        positive = float(row["positive_similarity"])
        negative = float(row["negative_similarity"])
        if not np.isfinite(positive) or not np.isfinite(negative):
            raise ValueError("non-finite relation similarity")
        contrast = positive - negative
        row["positive_similarity"] = positive
        row["negative_similarity"] = negative
        row["contrast"] = contrast
        row["F_credit"] = 1.0 if contrast > 0 else (0.0 if contrast == 0 else -1.0)
    return rows


def aggregate_relation_rows(rows, bootstrap_seed=BOOTSTRAP_SEED, bootstrap_draws=BOOTSTRAP_DRAWS):
    if not rows:
        raise ValueError("no relation rows")
    positive_groups = {}
    for row in rows:
        positive_groups.setdefault((row["anchor_id"], row["positive_relation_id"]), []).append(float(row["F_credit"]))
    positive_scores = {key: float(np.mean(values)) for key, values in positive_groups.items()}
    anchor_groups = {}
    for (anchor_id, _), value in positive_scores.items():
        anchor_groups.setdefault(anchor_id, []).append(value)
    anchor_scores = {anchor_id: float(np.mean(values)) for anchor_id, values in anchor_groups.items()}

    parent_of_anchor, unit_of_parent = {}, {}
    for row in rows:
        anchor_id, parent_id, unit_id = row["anchor_id"], row["source_parent_id"], row["source_dependency_unit_id"]
        if anchor_id in parent_of_anchor and parent_of_anchor[anchor_id] != parent_id:
            raise ValueError("anchor maps to multiple source parents")
        if parent_id in unit_of_parent and unit_of_parent[parent_id] != unit_id:
            raise ValueError("source parent maps to multiple dependency units")
        parent_of_anchor[anchor_id] = parent_id
        unit_of_parent[parent_id] = unit_id

    parent_groups = {}
    for anchor_id, value in anchor_scores.items():
        parent_groups.setdefault(parent_of_anchor[anchor_id], []).append(value)
    parent_scores = {parent_id: float(np.mean(values)) for parent_id, values in parent_groups.items()}
    unit_groups = {}
    for parent_id, value in parent_scores.items():
        unit_groups.setdefault(unit_of_parent[parent_id], []).append(value)
    unit_scores = {unit_id: float(np.mean(values)) for unit_id, values in unit_groups.items()}
    if not unit_scores:
        raise ValueError("no dependency-unit scores")

    score = float(np.mean(list(unit_scores.values())))
    rng = np.random.default_rng(bootstrap_seed)
    values = np.asarray(list(unit_scores.values()), dtype=np.float64)
    draws = np.mean(rng.choice(values, (bootstrap_draws, len(values)), replace=True), axis=1)
    return {
        "anchor_scores": anchor_scores,
        "parent_scores": parent_scores,
        "unit_scores": unit_scores,
        "parent_of_anchor": parent_of_anchor,
        "unit_of_parent": unit_of_parent,
        "structure_relation_group_F": score,
        "F": score,
        "F_95_CI": np.quantile(draws, [0.025, 0.975]).tolist(),
        "F_bootstrap": draws,
    }


def write_hierarchy_outputs(output_dir: Path, aggregation, prefix=""):
    stem = f"{prefix}_" if prefix else ""
    anchor_rows = [{
        "anchor_id": anchor_id,
        "source_parent_id": aggregation["parent_of_anchor"][anchor_id],
        "source_dependency_unit_id": aggregation["unit_of_parent"][aggregation["parent_of_anchor"][anchor_id]],
        "anchor_F": value,
    } for anchor_id, value in sorted(aggregation["anchor_scores"].items())]
    parent_rows = [{
        "source_parent_id": parent_id,
        "source_dependency_unit_id": aggregation["unit_of_parent"][parent_id],
        "source_parent_F": value,
    } for parent_id, value in sorted(aggregation["parent_scores"].items())]
    unit_rows = [{"source_dependency_unit_id": unit_id, "unit_F": value} for unit_id, value in sorted(aggregation["unit_scores"].items())]
    group_rows = [{"structure_relation_group_id": "RNA09_STRUCTURE_RELATION", "group_F": aggregation["structure_relation_group_F"]}]
    write_tsv(output_dir / f"{stem}INDIVIDUAL_SCORES.tsv", anchor_rows)
    write_tsv(output_dir / f"{stem}SOURCE_PARENT_SCORE_LEDGER.tsv", parent_rows)
    write_tsv(output_dir / f"{stem}SOURCE_UNIT_SCORE_LEDGER.tsv", unit_rows)
    write_tsv(output_dir / f"{stem}STRUCTURE_RELATION_GROUP_SCORE_LEDGER.tsv", group_rows)


def validate_refset(paths):
    fasta = read_fasta(paths["sequences"])
    object_rows = read_tsv(paths["object_map"])
    if not object_rows or not {"sequence_sha256", "object_id", "sequence_length_nt"}.issubset(object_rows[0]):
        raise ValueError("invalid endpoint object map schema")
    object_to_sequence = {}
    for row in object_rows:
        object_id, sequence_id = row["object_id"], row["sequence_sha256"]
        if object_id in object_to_sequence:
            raise ValueError(f"duplicate object_id: {object_id}")
        if sequence_id not in fasta:
            raise ValueError(f"object map sequence missing from FASTA: {sequence_id}")
        sequence = fasta[sequence_id]
        if hashlib.sha256(sequence.encode()).hexdigest() != sequence_id:
            raise ValueError(f"FASTA sequence identity mismatch: {sequence_id}")
        if len(sequence) != int(row["sequence_length_nt"]):
            raise ValueError(f"FASTA length mismatch: {object_id}")
        object_to_sequence[object_id] = sequence_id
    if set(fasta) != set(object_to_sequence.values()):
        raise ValueError("FASTA and object-map sequence universes differ")

    rows = read_tsv(paths["relations"])
    if not rows or len({row["contrast_id"] for row in rows}) != len(rows):
        raise ValueError("empty or duplicate contrast ledger")
    relation_ids, object_ids = {row["contrast_id"] for row in rows}, set(object_to_sequence)
    for row in rows:
        endpoints = {row[key] for key in ("anchor_id", "positive_partner_id", "negative_partner_id")}
        if len(endpoints) != 3 or not endpoints.issubset(object_ids):
            raise ValueError(f"invalid relation endpoints: {row['contrast_id']}")
        if not row["positive_relation_id"] or not row["negative_relation_id"]:
            raise ValueError(f"missing relation identity: {row['contrast_id']}")

    bundle_rows = read_tsv(paths["anchor_bundles"])
    if not bundle_rows or len({row["anchor_id"] for row in bundle_rows}) != len(bundle_rows):
        raise ValueError("empty or duplicate anchor bundle ledger")
    relation_by_anchor, bundled = {}, []
    for row in rows:
        relation_by_anchor.setdefault(row["anchor_id"], []).append(row)
    for bundle in bundle_rows:
        ids = [value for value in bundle["contrast_ids"].split(",") if value]
        bundled.extend(ids)
        if bundle["complete_bundle_status"] != "COMPLETE" or int(bundle["released_contrast_count"]) != len(ids):
            raise ValueError(f"incomplete anchor bundle: {bundle['anchor_id']}")
        expected = relation_by_anchor.get(bundle["anchor_id"], [])
        if set(ids) != {row["contrast_id"] for row in expected}:
            raise ValueError(f"anchor bundle contrast mismatch: {bundle['anchor_id']}")
        if any(row["source_parent_id"] != bundle["source_parent_id"] or row["source_dependency_unit_id"] != bundle["source_dependency_unit_id"] for row in expected):
            raise ValueError(f"anchor bundle hierarchy mismatch: {bundle['anchor_id']}")
    if len(bundled) != len(set(bundled)) or set(bundled) != relation_ids:
        raise ValueError("anchor bundles do not partition the contrast population")

    dependency_members = {}
    for row in read_tsv(paths["source_dependencies"]):
        unit_id = row["source_dependency_unit_id"]
        if unit_id in dependency_members:
            raise ValueError(f"duplicate dependency unit: {unit_id}")
        dependency_members[unit_id] = {value for value in row["member_object_ids"].split(",") if value}
    for row in rows:
        if row["source_dependency_unit_id"] not in dependency_members or row["source_parent_id"] not in dependency_members[row["source_dependency_unit_id"]]:
            raise ValueError(f"source parent/dependency mismatch: {row['contrast_id']}")

    declared_reuse = {(row["role"], row["object_id"]): row for row in read_tsv(paths["reuse"])}
    actual_reuse = {}
    for role, field in (("POSITIVE_TARGET", "positive_partner_id"), ("NEGATIVE_TARGET", "negative_partner_id")):
        for row in rows:
            key = (role, row[field])
            actual_reuse[key] = actual_reuse.get(key, 0) + 1
    if set(actual_reuse) != set(declared_reuse):
        raise ValueError("reuse ledger object universe mismatch")
    allocation_contracts = {row["allocation_contract_id"] for row in rows}
    for key, count in actual_reuse.items():
        declared = declared_reuse[key]
        if int(declared["released_reuse_count"]) != count or count > int(declared["frozen_cap"]):
            raise ValueError(f"reuse count mismatch: {key}")
        if declared["allocation_contract_id"] not in allocation_contracts:
            raise ValueError(f"reuse allocation contract mismatch: {key}")

    manifest = json.loads(paths["population_manifest"].read_text(encoding="utf-8"))
    actual_counts = {
        "anchor_bundles": len(bundle_rows), "anchors": len({row["anchor_id"] for row in rows}),
        "contrasts": len(rows), "negative_targets": len({row["negative_partner_id"] for row in rows}),
        "positive_targets": len({row["positive_partner_id"] for row in rows}),
        "source_dependency_units": len({row["source_dependency_unit_id"] for row in rows}),
        "source_parents": len({row["source_parent_id"] for row in rows}),
    }
    if actual_counts != manifest["planned_counts"]:
        raise ValueError(f"planned population counts mismatch: {actual_counts}")
    if allocation_contracts != {manifest["allocation_contract_id"]}:
        raise ValueError("population allocation contract mismatch")
    return object_to_sequence, rows, manifest, actual_counts


def main():
    parser = argparse.ArgumentParser(description="Portable RNA09 formal scorer")
    root = Path(__file__).resolve().parents[1]
    defaults = {
        "sequences": root / "refset/scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta",
        "object_map": root / "refset/scoring_input_v01/RNA09_ENDPOINT_OBJECT_MAP.tsv",
        "relations": root / "refset/FORMAL_CONTRAST_LEDGER.tsv",
        "anchor_bundles": root / "refset/ANCHOR_BUNDLE_LEDGER.tsv",
        "source_dependencies": root / "refset/SOURCE_DEPENDENCY_LEDGER.tsv",
        "reuse": root / "refset/REUSE_LEDGER.tsv",
        "population_manifest": root / "refset/PLANNED_POPULATION_MANIFEST.json",
    }
    parser.add_argument("--embeddings", type=Path, required=True, help="NPZ with sequence_sha256 and vectors")
    for key, default in defaults.items():
        parser.add_argument(f"--{key.replace('_', '-')}", type=Path, default=default)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.output_dir.exists():
        raise FileExistsError(f"output directory already exists: {args.output_dir}")
    data = np.load(args.embeddings, allow_pickle=False)
    if set(data.files) != {"sequence_sha256", "vectors"}:
        raise ValueError("embedding NPZ must contain only sequence_sha256 and vectors")
    sequence_ids = [str(value) for value in data["sequence_sha256"]]
    raw_vectors = data["vectors"]
    if np.issubdtype(raw_vectors.dtype, np.bool_) or not np.issubdtype(raw_vectors.dtype, np.number) or np.issubdtype(raw_vectors.dtype, np.complexfloating):
        raise ValueError("vectors must use a real numeric dtype")
    vectors = np.asarray(raw_vectors, dtype=np.float64)
    if len(sequence_ids) != len(set(sequence_ids)) or vectors.ndim != 2 or vectors.shape[1] <= 0 or vectors.shape[0] != len(sequence_ids):
        raise ValueError("invalid embedding identity or dimensions")
    if not np.isfinite(vectors).all():
        raise ValueError("non-finite embeddings")
    vector_by_sequence = dict(zip(sequence_ids, vectors))

    paths = {key: Path(getattr(args, key)) for key in defaults}
    object_to_sequence, rows, population_manifest, population_counts = validate_refset(paths)
    if set(sequence_ids) != set(object_to_sequence.values()):
        raise ValueError("embedding IDs must exactly cover the released endpoint sequence set")
    for row in rows:
        anchor, positive, negative = (vector_by_sequence[object_to_sequence[row[key]]] for key in ("anchor_id", "positive_partner_id", "negative_partner_id"))
        row["positive_similarity"], row["negative_similarity"] = cosine(anchor, positive), cosine(anchor, negative)
    assign_credits(rows)
    aggregation = aggregate_relation_rows(rows)

    args.output_dir.mkdir(parents=True)
    leaf_path = args.output_dir / "LEAF_LEDGER.tsv"
    write_tsv(leaf_path, rows)
    write_hierarchy_outputs(args.output_dir, aggregation)
    np.save(args.output_dir / "F_BOOTSTRAP_MEANS.npy", aggregation["F_bootstrap"])
    actual_hashes = {key: digest(path) for key, path in paths.items()}
    registered_hashes = {key: digest(path) for key, path in defaults.items()}
    result = {
        "status": "PASS", "constraint_id": "RNA09", "release": "v0.1.3_F_ONLY", "formal_result": True,
        "registered_refset_match": actual_hashes == registered_hashes, "refset_input_sha256": actual_hashes,
        "planned_population_identity_sha256": population_manifest["planned_population_identity_sha256"],
        "records": len(vector_by_sequence), "contrasts": len(rows), "anchors": population_counts["anchors"],
        "source_parents": population_counts["source_parents"], "source_dependency_units": len(aggregation["unit_scores"]),
        "structure_relation_groups": 1, "F": aggregation["F"], "F_95_CI": aggregation["F_95_CI"],
        "formal_interface": "F_ONLY", "score_transform_version": "RNA09_F_ONLY_v1", "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "leaf_ledger_sha256": digest(leaf_path),
    }
    (args.output_dir / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    outputs = [path for path in args.output_dir.iterdir() if path.is_file()]
    manifest = {"files": [{"path": path.name, "bytes": path.stat().st_size, "sha256": digest(path)} for path in sorted(outputs)]}
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
