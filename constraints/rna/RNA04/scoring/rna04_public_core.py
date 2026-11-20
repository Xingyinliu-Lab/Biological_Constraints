from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


STANDARD_CODE = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*", "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records[current] = "".join(chunks).upper()
            current = line[1:].split()[0]
            if not current or current in records:
                raise ValueError("FASTA IDs must be nonempty and unique")
            chunks = []
        elif current is None:
            raise ValueError("sequence before first FASTA header")
        else:
            chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks).upper()
    if not records:
        raise ValueError("empty endpoint FASTA")
    return records


def validate_refset(rows: list[dict], fasta_path: Path) -> dict:
    sequences = read_fasta(fasta_path)
    required = {row[key] for row in rows for key in ("wt_sequence_sha256", "alt_sequence_sha256")}
    if set(sequences) != required:
        raise ValueError(f"FASTA coverage mismatch; missing={len(required-set(sequences))}, extra={len(set(sequences)-required)}")
    parent_cluster: dict[str, str] = {}
    parent_wt_sequence: dict[str, str] = {}
    site_identity: dict[str, tuple[str, str, int, str]] = {}
    seen_site_variants: set[tuple[str, str]] = set()
    semantic_rows: list[dict] = []
    for sequence_id, sequence in sequences.items():
        if not sequence or set(sequence) - set("ACGU") or len(sequence) % 3:
            raise ValueError(f"invalid complete coding RNA sequence: {sequence_id}")
        if hashlib.sha256(sequence.encode()).hexdigest() != sequence_id:
            raise ValueError(f"FASTA sequence identity mismatch: {sequence_id}")
    for index, row in enumerate(rows):
        wt_id, alt_id = row["wt_sequence_sha256"], row["alt_sequence_sha256"]
        wt, alt = sequences[wt_id], sequences[alt_id]
        position = int(row["position0"])
        if not 0 <= position < len(wt) or len(wt) != len(alt):
            raise ValueError(f"invalid WT/ALT coordinates at row {index}")
        differences = [offset for offset, pair in enumerate(zip(wt, alt)) if pair[0] != pair[1]]
        if differences != [position] or row.get("hamming") != "1":
            raise ValueError(f"row {index} is not the declared single-nucleotide variant")
        if alt != row["alt_sequence"] or row["alt"] != alt[position]:
            raise ValueError(f"ALT sequence annotation mismatch at row {index}")
        start = position - position % 3
        wt_codon, alt_codon = wt[start:start+3], alt[start:start+3]
        if wt_codon != row["wt_codon"].replace("T", "U") or alt_codon != row["alt_codon"].replace("T", "U"):
            raise ValueError(f"codon annotation mismatch at row {index}")
        wt_aa, alt_aa = STANDARD_CODE[wt_codon], STANDARD_CODE[alt_codon]
        expected_class = "synonymous" if wt_aa == alt_aa else "missense"
        if "*" in (wt_aa, alt_aa) or row["wt_amino_acid"] != wt_aa or row["alt_amino_acid"] != alt_aa or row["mutation_class"] != expected_class:
            raise ValueError(f"Codon consequence annotation mismatch at row {index}")
        parent, cluster, site = row["parent_object_id"], row["outer_cluster_id"], row["site_id"]
        site_variant = (site, alt_id)
        if site_variant in seen_site_variants:
            raise ValueError("duplicate ALT endpoint within site")
        seen_site_variants.add(site_variant)
        if parent in parent_cluster and parent_cluster[parent] != cluster:
            raise ValueError("parent maps to multiple dependency clusters")
        parent_cluster[parent] = cluster
        if parent in parent_wt_sequence and parent_wt_sequence[parent] != wt_id:
            raise ValueError("parent maps to multiple WT sequence identities")
        parent_wt_sequence[parent] = wt_id
        identity = (parent, cluster, position, wt_id)
        if site in site_identity and site_identity[site] != identity:
            raise ValueError("site maps to inconsistent biological identities")
        site_identity[site] = identity
        semantic_rows.append({key: row[key] for key in (
            "site_id", "parent_object_id", "outer_cluster_id", "position0", "wt_sequence_sha256",
            "alt_sequence_sha256", "mutation_class", "wt_codon", "alt_codon",
        )})
    identity_payload = {
        "sequences": sorted(sequences),
        "rows": sorted(semantic_rows, key=lambda value: tuple(value.values())),
    }
    identity = hashlib.sha256(json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"sequences": sequences, "population_identity_sha256": identity}


def pair_F_credit(left: float, right: float) -> float:
    """Formal F-only credit for one missense-versus-synonymous comparison."""
    return 1.0 if left > right else -1.0 if left < right else 0.0


def site_scores(rows: list[dict], scores: np.ndarray) -> tuple[dict[str, float], list[dict], list[dict]]:
    grouped: dict[str, dict[str, list[tuple[dict, float]]]] = defaultdict(lambda: {"missense": [], "synonymous": []})
    identity: dict[str, tuple[str, str]] = {}
    for row, score in zip(rows, scores):
        grouped[row["site_id"]][row["mutation_class"]].append((row, float(score)))
        current = (row["parent_object_id"], row["outer_cluster_id"])
        if row["site_id"] in identity and identity[row["site_id"]] != current:
            raise ValueError("site identity maps to multiple parent/cluster units")
        identity[row["site_id"]] = current
    values: dict[str, float] = {}
    records: list[dict] = []
    pair_records: list[dict] = []
    for site_id in sorted(grouped):
        classes = grouped[site_id]
        if not classes["missense"] or not classes["synonymous"]:
            raise ValueError(f"incomplete site {site_id}")
        pairs = [
            (missense_row, missense_score, synonymous_row, synonymous_score)
            for missense_row, missense_score in classes["missense"]
            for synonymous_row, synonymous_score in classes["synonymous"]
        ]
        credits = [pair_F_credit(missense_score, synonymous_score) for _, missense_score, _, synonymous_score in pairs]
        site_F = float(np.mean(credits))
        values[site_id] = site_F
        parent_id, cluster_id = identity[site_id]
        records.append({
            "site_id": site_id, "parent_object_id": parent_id, "outer_cluster_id": cluster_id,
            "missense_count": len(classes["missense"]), "synonymous_count": len(classes["synonymous"]),
            "pairwise_comparisons": len(credits), "individual_F": format(site_F, ".17g"),
        })
        for missense_row, missense_score, synonymous_row, synonymous_score in pairs:
            margin = missense_score - synonymous_score
            pair_id = hashlib.sha256(
                f"{site_id}|{missense_row['alt_sequence_sha256']}|{synonymous_row['alt_sequence_sha256']}".encode()
            ).hexdigest()
            pair_records.append({
                "pair_id": pair_id,
                "site_id": site_id,
                "parent_object_id": parent_id,
                "outer_cluster_id": cluster_id,
                "missense_alt_sequence_sha256": missense_row["alt_sequence_sha256"],
                "missense_score": format(missense_score, ".17g"),
                "synonymous_alt_sequence_sha256": synonymous_row["alt_sequence_sha256"],
                "synonymous_score": format(synonymous_score, ".17g"),
                "margin": format(margin, ".17g"),
                "tie": str(margin == 0.0).lower(),
                "F_credit": format(pair_F_credit(missense_score, synonymous_score), ".17g"),
            })
    return values, records, pair_records


def aggregate(rows: list[dict], site_values: dict[str, float]) -> tuple[float, dict[str, float], dict[str, float]]:
    parent_sites: dict[str, list[float]] = defaultdict(list)
    parent_cluster: dict[str, str] = {}
    seen_parent_site: set[tuple[str, str]] = set()
    for row in rows:
        site_id = row["site_id"]
        parent_id = row["parent_object_id"]
        if site_id not in site_values or (parent_id, site_id) in seen_parent_site:
            continue
        seen_parent_site.add((parent_id, site_id))
        parent_sites[parent_id].append(site_values[site_id])
        parent_cluster[parent_id] = row["outer_cluster_id"]
    parent_values = {parent: float(np.mean(values)) for parent, values in parent_sites.items()}
    cluster_parents: dict[str, list[float]] = defaultdict(list)
    for parent, value in parent_values.items():
        cluster_parents[parent_cluster[parent]].append(value)
    cluster_values = {cluster: float(np.mean(values)) for cluster, values in cluster_parents.items()}
    return float(np.mean(list(cluster_values.values()))), cluster_values, parent_values
