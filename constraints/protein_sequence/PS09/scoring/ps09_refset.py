"""Semantic closure validation for the released PS09 RefSet."""

from __future__ import annotations

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_MEMBERSHIP_FIELDS = {
    "block_id",
    "homology_cluster_id",
    "uniprot_accession",
    "sequence_length",
    "uniprot_position_1based",
    "amino_acid",
    "relative_position_bin",
    "state",
}
FORMAL_STATES = {"DISORDER", "ORDER_COORDINATE_SUPPORTED"}


def endpoint_id(row: dict[str, str]) -> str:
    return f"{row['block_id']}|{row['state']}|{row['uniprot_position_1based']}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"missing TSV header: {path}")
        return list(reader)


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    identity: str | None = None
    residues: list[str] = []
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identity is not None:
                    if identity in sequences:
                        raise ValueError(f"duplicate FASTA accession: {identity}")
                    sequences[identity] = "".join(residues)
                identity = line[1:].split()[0]
                if not identity:
                    raise ValueError("empty FASTA accession")
                residues = []
            else:
                if identity is None:
                    raise ValueError("FASTA sequence precedes first header")
                residues.append(line.upper())
    if identity is not None:
        if identity in sequences:
            raise ValueError(f"duplicate FASTA accession: {identity}")
        sequences[identity] = "".join(residues)
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError("empty FASTA or FASTA sequence")
    return sequences


def validate_refset(
    *, membership: Path, blocks: Path, clusters: Path, fasta: Path, folds: Path
) -> dict[str, object]:
    membership_rows = read_tsv(membership)
    if not membership_rows:
        raise ValueError("empty membership table")
    if not REQUIRED_MEMBERSHIP_FIELDS.issubset(membership_rows[0]):
        missing = sorted(REQUIRED_MEMBERSHIP_FIELDS - set(membership_rows[0]))
        raise ValueError(f"membership fields missing: {missing}")
    endpoint_ids = [endpoint_id(row) for row in membership_rows]
    if len(endpoint_ids) != len(set(endpoint_ids)):
        raise ValueError("duplicate formal endpoint ID")

    sequences = read_fasta(fasta)
    membership_accessions = {row["uniprot_accession"] for row in membership_rows}
    if membership_accessions != set(sequences):
        raise ValueError("membership/FASTA accession closure failure")

    block_members: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in membership_rows:
        accession = row["uniprot_accession"]
        sequence = sequences[accession]
        try:
            declared_length = int(row["sequence_length"])
            position = int(row["uniprot_position_1based"])
        except ValueError as exc:
            raise ValueError("non-integer sequence length or residue coordinate") from exc
        if declared_length != len(sequence):
            raise ValueError(f"sequence length mismatch: {accession}")
        if not 1 <= position <= len(sequence):
            raise ValueError(f"residue coordinate out of range: {endpoint_id(row)}")
        if row["amino_acid"] != sequence[position - 1]:
            raise ValueError(f"residue identity mismatch: {endpoint_id(row)}")
        if row["state"] not in FORMAL_STATES:
            raise ValueError(f"invalid endpoint state: {row['state']}")
        block_members[row["block_id"]].append(row)

    block_rows = read_tsv(blocks)
    required_block_fields = {
        "block_id",
        "uniprot_accession",
        "amino_acid",
        "relative_position_bin",
        "disorder_endpoint_count",
        "order_endpoint_count",
    }
    if not block_rows or not required_block_fields.issubset(block_rows[0]):
        raise ValueError("invalid matched-block table")
    block_by_id = {row["block_id"]: row for row in block_rows}
    if len(block_by_id) != len(block_rows) or set(block_by_id) != set(block_members):
        raise ValueError("membership/matched-block closure failure")
    block_clusters: dict[str, str] = {}
    for block_id, rows in block_members.items():
        attributes = {
            (
                row["uniprot_accession"],
                row["amino_acid"],
                row["relative_position_bin"],
                row["homology_cluster_id"],
            )
            for row in rows
        }
        if len(attributes) != 1:
            raise ValueError(f"inconsistent biological identity within block: {block_id}")
        accession, amino_acid, relative_bin, cluster_id = next(iter(attributes))
        block_clusters[block_id] = cluster_id
        inventory = block_by_id[block_id]
        counts = Counter(row["state"] for row in rows)
        expected = (
            inventory["uniprot_accession"],
            inventory["amino_acid"],
            inventory["relative_position_bin"],
            int(inventory["disorder_endpoint_count"]),
            int(inventory["order_endpoint_count"]),
        )
        observed = (
            accession,
            amino_acid,
            relative_bin,
            counts["DISORDER"],
            counts["ORDER_COORDINATE_SUPPORTED"],
        )
        if observed != expected:
            raise ValueError(f"matched-block membership mismatch: {block_id}")

    cluster_pairs: dict[str, str] = {}
    with clusters.open(encoding="utf-8-sig", newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) != 2 or not all(fields):
                raise ValueError("MMseqs2 cluster inventory must contain two columns")
            representative, accession = fields
            if accession in cluster_pairs:
                raise ValueError(f"duplicate cluster member: {accession}")
            cluster_pairs[accession] = representative
    if set(cluster_pairs) != set(sequences):
        raise ValueError("cluster inventory/FASTA accession closure failure")
    for row in membership_rows:
        expected_cluster = f"mmseqs2:{cluster_pairs[row['uniprot_accession']]}"
        if row["homology_cluster_id"] != expected_cluster:
            raise ValueError(f"homology cluster mismatch: {endpoint_id(row)}")

    fold_rows = read_tsv(folds)
    required_fold_fields = {"homology_cluster_id", "fold", "block_count"}
    if not fold_rows or not required_fold_fields.issubset(fold_rows[0]):
        raise ValueError("invalid fold-assignment table")
    fold_map = {row["homology_cluster_id"]: int(row["fold"]) for row in fold_rows}
    membership_clusters = set(block_clusters.values())
    if len(fold_map) != len(fold_rows) or set(fold_map) != membership_clusters:
        raise ValueError("membership/fold cluster closure failure")
    if set(fold_map.values()) != set(range(5)):
        raise ValueError("invalid frozen five-fold assignment")
    observed_block_counts = Counter(block_clusters.values())
    for row in fold_rows:
        if int(row["block_count"]) != observed_block_counts[row["homology_cluster_id"]]:
            raise ValueError(f"fold block-count mismatch: {row['homology_cluster_id']}")

    return {
        "membership_rows": membership_rows,
        "folds": fold_map,
        "counts": {
            "proteins": len(sequences),
            "endpoints": len(membership_rows),
            "blocks": len(block_members),
            "homology_clusters": len(membership_clusters),
        },
    }
