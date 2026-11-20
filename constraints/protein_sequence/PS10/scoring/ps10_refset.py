"""Semantic validation for PS10 records and study-aware fold assignments."""

from __future__ import annotations

import csv
import hashlib
import re
from collections import Counter
from pathlib import Path


FORMAL_CLASSES = {"AMYLOID", "NON_AMYLOID"}
FORMAL_ROLES = {"train", "test", "purged"}
PEPTIDE_PATTERN = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]+")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"missing TSV header: {path}")
        return list(reader)


def study_units(row: dict[str, str]) -> set[str]:
    units = {value for value in row["source_study_units"].split(";") if value}
    if not units:
        raise ValueError(f"record has no source study: {row['record_id']}")
    return units


def validate_records(records_path: Path) -> list[dict[str, str]]:
    rows = read_tsv(records_path)
    required = {
        "record_id",
        "normalized_peptide_identity",
        "peptide_identity_sha256",
        "experimental_classification",
        "length",
        "source_study_units",
        "dependency_component_id",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("invalid PS10 record-table schema")
    ids = [row["record_id"] for row in rows]
    sequences = [row["normalized_peptide_identity"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate frozen record ID")
    if len(sequences) != len(set(sequences)):
        raise ValueError("duplicate normalized peptide identity")
    for row in rows:
        sequence = row["normalized_peptide_identity"]
        if not PEPTIDE_PATTERN.fullmatch(sequence):
            raise ValueError(f"invalid normalized peptide sequence: {row['record_id']}")
        if int(row["length"]) != len(sequence):
            raise ValueError(f"peptide length mismatch: {row['record_id']}")
        if hashlib.sha256(sequence.encode()).hexdigest() != row["peptide_identity_sha256"]:
            raise ValueError(f"peptide identity hash mismatch: {row['record_id']}")
        if row["experimental_classification"] not in FORMAL_CLASSES:
            raise ValueError(
                f"invalid experimental classification: {row['experimental_classification']}"
            )
        if not row["dependency_component_id"]:
            raise ValueError(f"missing dependency component: {row['record_id']}")
        study_units(row)
    return rows


def validate_refset(records_path: Path, assignments_path: Path) -> dict[str, object]:
    records = validate_records(records_path)
    assignments = read_tsv(assignments_path)
    required = {"fold_id", "record_id", "role"}
    if not assignments or not required.issubset(assignments[0]):
        raise ValueError("invalid PS10 assignment schema")
    by_id = {row["record_id"]: row for row in records}
    folds = sorted({row["fold_id"] for row in assignments})
    expected_folds = [f"STUDY-{index:02d}" for index in range(1, 23)]
    if folds != expected_folds:
        raise ValueError("PS10 requires the frozen 22 study folds")
    pairs = [(row["fold_id"], row["record_id"]) for row in assignments]
    if len(pairs) != len(set(pairs)) or len(assignments) != len(records) * 22:
        raise ValueError("assignment table is not one role per record per fold")
    if {row["record_id"] for row in assignments} != set(by_id):
        raise ValueError("assignment/record population closure failure")
    if any(row["role"] not in FORMAL_ROLES for row in assignments):
        raise ValueError("invalid fold role")

    test_occurrences: Counter[str] = Counter()
    heldout_studies: dict[str, str] = {}
    for fold_id in folds:
        rows = [row for row in assignments if row["fold_id"] == fold_id]
        train = [row["record_id"] for row in rows if row["role"] == "train"]
        test = [row["record_id"] for row in rows if row["role"] == "test"]
        if not train or not test:
            raise ValueError(f"empty train or test population: {fold_id}")
        train_classes = {by_id[identity]["experimental_classification"] for identity in train}
        test_classes = {by_id[identity]["experimental_classification"] for identity in test}
        if train_classes != FORMAL_CLASSES or test_classes != FORMAL_CLASSES:
            raise ValueError(f"fold lacks both formal classes: {fold_id}")
        common_studies = set.intersection(*(study_units(by_id[identity]) for identity in test))
        if len(common_studies) != 1:
            raise ValueError(f"held-out study is not uniquely defined: {fold_id}")
        heldout_study = next(iter(common_studies))
        if any(heldout_study in study_units(by_id[identity]) for identity in train):
            raise ValueError(f"held-out study leaked into training: {fold_id}")
        heldout_studies[fold_id] = heldout_study
        test_occurrences.update(test)
    if any(count != 1 for count in test_occurrences.values()):
        raise ValueError("a record is held out in more than one primary fold")

    return {
        "records": records,
        "assignments": assignments,
        "counts": {
            "records": len(records),
            "folds": len(folds),
            "assignment_rows": len(assignments),
            "heldout_records": len(test_occurrences),
        },
        "heldout_studies": heldout_studies,
    }
