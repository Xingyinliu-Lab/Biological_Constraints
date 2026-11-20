#!/usr/bin/env python3
"""PS08 population-schema validation and provenance helpers."""
from __future__ import annotations

import hashlib
import gzip
import json
import re
from pathlib import Path


REGISTERED_POPULATION_FILE_SHA256 = "798373b003ba4af39fa13f8c292351154673272100d83d8102a07d83e9a4ebe9"
REGISTERED_POPULATION_ID = "2ee0bbde63d2ed1bcec5b0c12c00b8a10385cc98523a9aa8966ae1a4fd5a9fb3"
ID_FIELDS = (
    "leaf_id",
    "first_edit_id",
    "source_parent_id",
    "source_unit_id",
    "x0_object_id",
    "x1_object_id",
    "x2_object_id",
)
SEQUENCE_HASH_FIELDS = ("x0_sequence_sha256", "x1_sequence_sha256", "x2_sequence_sha256")
REQUIRED_FIELDS = (*ID_FIELDS, "nested_extension_id", *SEQUENCE_HASH_FIELDS)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    sequences: dict[str, str] = {}
    object_id: str | None = None
    chunks: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if object_id is not None:
                    if object_id in sequences:
                        raise ValueError(f"duplicate FASTA object ID: {object_id}")
                    sequences[object_id] = "".join(chunks)
                object_id = line[1:].split()[0]
                chunks = []
            else:
                if object_id is None:
                    raise ValueError("FASTA sequence occurs before the first header")
                chunks.append(line.upper())
    if object_id is not None:
        if object_id in sequences:
            raise ValueError(f"duplicate FASTA object ID: {object_id}")
        sequences[object_id] = "".join(chunks)
    if not sequences:
        raise ValueError("endpoint FASTA is empty")
    for key, sequence in sequences.items():
        if not sequence or set(sequence) - AMINO_ACIDS:
            raise ValueError(f"invalid canonical amino-acid sequence for {key}")
    return sequences


def validate_sequence_relations(members: list[dict], fasta_path: Path) -> dict:
    sequences = load_fasta(fasta_path)
    expected_ids = {
        row[field]
        for row in members
        for field in ("x0_object_id", "x1_object_id", "x2_object_id")
    }
    if set(sequences) != expected_ids:
        raise ValueError("endpoint FASTA object universe differs from population")
    for index, row in enumerate(members):
        ladder = [sequences[row[f"x{i}_object_id"]] for i in range(3)]
        for i, sequence in enumerate(ladder):
            digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            if digest != row[f"x{i}_sequence_sha256"]:
                raise ValueError(f"population member {index} sequence SHA256 mismatch")
        if len({len(sequence) for sequence in ladder}) != 1:
            raise ValueError(f"population member {index} endpoints have different lengths")
        h01 = sum(a != b for a, b in zip(ladder[0], ladder[1]))
        h12 = sum(a != b for a, b in zip(ladder[1], ladder[2]))
        h02 = sum(a != b for a, b in zip(ladder[0], ladder[2]))
        if (h01, h12, h02) != (1, 1, 2):
            raise ValueError(f"population member {index} violates the PS08 Hamming 1/1/2 relation")
    return {
        "endpoint_fasta_sha256": sha256(fasta_path),
        "endpoint_count": len(sequences),
        "hamming_relation": "H(x0,x1)=1;H(x1,x2)=1;H(x0,x2)=2",
        "status": "PASS",
    }


def load_population(path: Path, fasta_path: Path | None = None) -> tuple[list[dict], dict]:
    population = json.loads(path.read_text(encoding="utf-8-sig"))
    if population.get("kind") != "PS08_FROZEN_POPULATION_MANIFEST":
        raise ValueError("PS08 population kind is invalid")
    members = population.get("canonical_members")
    if not isinstance(members, list) or not members:
        raise ValueError("PS08 population must contain canonical_members")
    for index, row in enumerate(members):
        if not isinstance(row, dict) or any(not isinstance(row.get(field), str) or not row[field] for field in REQUIRED_FIELDS):
            raise ValueError(f"population member {index} lacks a required non-empty string field")
        if any(not SHA256_RE.fullmatch(row[field]) for field in SEQUENCE_HASH_FIELDS):
            raise ValueError(f"population member {index} has an invalid sequence SHA256")
        if len({row["x0_object_id"], row["x1_object_id"], row["x2_object_id"]}) != 3:
            raise ValueError(f"population member {index} does not contain three distinct endpoints")
    for field in ("leaf_id", "first_edit_id", "nested_extension_id"):
        values = [row[field] for row in members]
        if len(values) != len(set(values)):
            raise ValueError(f"duplicate {field} in PS08 population")
    object_ids = [row[field] for row in members for field in ("x0_object_id", "x1_object_id", "x2_object_id")]
    if len(object_ids) != len(set(object_ids)):
        raise ValueError("PS08 endpoint object IDs are not globally unique")
    parent_to_unit: dict[str, str] = {}
    first_to_parent: dict[str, tuple[str, str]] = {}
    for row in members:
        parent, unit = row["source_parent_id"], row["source_unit_id"]
        if parent in parent_to_unit and parent_to_unit[parent] != unit:
            raise ValueError("one source parent maps to multiple source units")
        parent_to_unit[parent] = unit
        first = row["first_edit_id"]
        mapping = (parent, unit)
        if first in first_to_parent and first_to_parent[first] != mapping:
            raise ValueError("one first edit maps to inconsistent parent or unit identities")
        first_to_parent[first] = mapping
    semantic_records = [{field: row[field] for field in REQUIRED_FIELDS} for row in members]
    semantic_records.sort(key=lambda row: tuple(row[field] for field in REQUIRED_FIELDS))
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic_records, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()
    actual_file_sha256 = sha256(path)
    declared_identity = population.get("population_identity_sha256")
    if not isinstance(declared_identity, str) or not SHA256_RE.fullmatch(declared_identity):
        raise ValueError("population_identity_sha256 must be a lowercase SHA256")
    if declared_identity != semantic_sha256:
        raise ValueError("declared population identity does not match canonical_members")
    provenance = {
        "population_file_sha256": actual_file_sha256,
        "registered_population_file_sha256": REGISTERED_POPULATION_FILE_SHA256,
        "registered_population_bytes_match": actual_file_sha256 == REGISTERED_POPULATION_FILE_SHA256,
        "declared_population_identity_sha256": declared_identity,
        "registered_population_identity_sha256": REGISTERED_POPULATION_ID,
        "declared_population_identity_match": True,
        "registered_population_identity_match": semantic_sha256 == REGISTERED_POPULATION_ID,
        "semantic_population_sha256": semantic_sha256,
        "leaf_count": len(members),
        "object_count": len(object_ids),
    }
    if fasta_path is not None:
        provenance["sequence_validation"] = validate_sequence_relations(members, fasta_path)
    return members, provenance
