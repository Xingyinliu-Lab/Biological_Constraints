"""Portable form of the executed RNA05 ViennaRNA comparator adapter.

The measurement and aggregation logic is preserved from the original executed
2026-10-04 script. Only author-machine paths were replaced by CLI arguments.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path


CODE_DIR = Path(__file__).resolve().parent
COMPARATOR_DIR = CODE_DIR.parent
CONSTRAINT_DIR = COMPARATOR_DIR.parents[1]
sys.path.insert(0, str(COMPARATOR_DIR))

from vienna_partition_pairing_propensity import pairing_propensity  # noqa: E402
import frozen_aggregation_core as frozen  # noqa: E402


def canon(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha_obj(value: object) -> str:
    return hashlib.sha256(canon(value).encode()).hexdigest()


def sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                sequence = "".join(chunks)
                sequences[sequence_sha256(sequence)] = sequence
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip().upper())
    if name is not None:
        sequence = "".join(chunks)
        sequences[sequence_sha256(sequence)] = sequence
    return sequences


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta",
        type=Path,
        default=CONSTRAINT_DIR / "refset" / "RNA05_06_unique_inputs.fasta",
    )
    parser.add_argument(
        "--endpoint-ledger",
        type=Path,
        default=CONSTRAINT_DIR / "refset" / "FORMAL_8BIN_ENDPOINT_LEDGER.tsv",
    )
    parser.add_argument(
        "--pair-maps",
        type=Path,
        default=COMPARATOR_DIR / "method" / "complete_pair_state_maps.tsv",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    endpoints = read_tsv(args.endpoint_ledger)
    maps: dict[tuple[str, str, str], set[int]] = {}
    for row in read_tsv(args.pair_maps):
        pairs = json.loads(row["v1_eligible_positive_pairs_json"])
        positions = {int(value) for pair in pairs for value in pair[:2]}
        maps[(row["record_id"], row["structure_map_id"], row["condition_id"])] = positions
    sequences = read_fasta(args.fasta)

    block_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in endpoints:
        block_rows[row["block_id"]].append(row)
    qcache: dict[str, list[float]] = {}
    sequence_audit: list[dict[str, object]] = []
    for sequence_hash in sorted({row["sequence_sha256"] for row in endpoints}):
        sequence = sequences.get(sequence_hash)
        if sequence is None:
            raise ValueError(f"missing raw sequence for {sequence_hash}")
        propensity = pairing_propensity(sequence)
        if len(propensity) != len(sequence):
            raise ValueError(f"Vienna coordinate mismatch {sequence_hash}")
        qcache[sequence_hash] = propensity
        sequence_audit.append(
            {
                "sequence_sha256": sequence_hash,
                "length_nt": len(sequence),
                "pairing_propensity_sha256": sha_obj(propensity),
            }
        )

    records: list[dict[str, object]] = []
    compatible: list[str] = []
    out_of_scope: list[dict[str, object]] = []
    for block_id, rows in sorted(block_rows.items()):
        key = (
            rows[0]["record_id"],
            rows[0]["primary_structure_map_id"],
            rows[0]["primary_condition_id"],
        )
        support = maps.get(key)
        if support is None:
            raise ValueError(f"missing frozen map for {key!r}")
        paired = [row for row in rows if row["paired_state"] == "paired"]
        unpaired = [row for row in rows if row["paired_state"] == "unpaired"]
        if not paired or not unpaired:
            raise ValueError(f"incomplete frozen block {block_id}")
        unsupported = sorted(
            {
                int(row["position0"])
                for row in paired
                if int(row["position0"]) not in support
            }
        )
        if unsupported:
            out_of_scope.append(
                {
                    "block_id": block_id,
                    "reason": "PAIRED_TRUTH_OUTSIDE_VIENNA_SUPPORTED_SECONDARY_PAIR_UNIVERSE",
                    "unsupported_paired_position0": unsupported,
                }
            )
            continue
        compatible.append(block_id)
        paired_scores = [qcache[row["sequence_sha256"]][int(row["position0"])] for row in paired]
        unpaired_scores = [qcache[row["sequence_sha256"]][int(row["position0"])] for row in unpaired]
        mean_paired_score = math.fsum(paired_scores) / len(paired_scores)
        mean_unpaired_score = math.fsum(unpaired_scores) / len(unpaired_scores)
        credit = (
            1.0
            if mean_paired_score > mean_unpaired_score
            else -1.0 if mean_paired_score < mean_unpaired_score else 0.0
        )
        leaf_id = sha_obj(
            {
                "block_id": block_id,
                "credit_definition": "mean paired q versus mean explicitly-unpaired q",
            }
        )
        member = {
            "leaf_id": leaf_id,
            "class_contrast_id": leaf_id,
            "matched_block_id": block_id,
            "parent_id": paired[0]["parent_id"],
            "dependency_cluster_id": paired[0]["outer_cluster_id"],
        }
        records.append(
            {
                "leaf_id": leaf_id,
                "member": member,
                "F_credit": credit,
                "paired_position_count": len(paired),
                "unpaired_position_count": len(unpaired),
                "mean_paired_score": mean_paired_score,
                "mean_unpaired_score": mean_unpaired_score,
            }
        )
    if not records:
        raise ValueError("no comparator-compatible complete blocks")

    members = [record["member"] for record in records]
    levels = ["class_contrast_id", "matched_block_id", "parent_id", "dependency_cluster_id"]
    frozen.validate_members(members, levels)
    population_identity = sha_obj(sorted(members, key=canon))
    score, trace = frozen.aggregate(
        records,
        levels,
        ["CLASS_CONTRAST", "MATCHED_BLOCK", "PARENT", "DEPENDENCY_CLUSTER", "OVERALL"],
        population_identity,
    )
    cluster_values = [node["value"] for node in trace if node["level"] == "DEPENDENCY_CLUSTER"]
    rng = random.Random(20261004)
    draws = sorted(
        math.fsum(rng.choice(cluster_values) for _ in cluster_values) / len(cluster_values)
        for _ in range(20000)
    )
    confidence_interval = [
        draws[int(0.025 * (len(draws) - 1))],
        draws[int(0.975 * (len(draws) - 1))],
    ]
    compatibility = {
        "comparator_id": "CMP-RNA05-VIENNA-PARTITION-PROPENSITY",
        "input_boundary": "raw RNA sequence only; no map, coordinates, labels, constraints, or chemical reactivity supplied to ViennaRNA",
        "compatibility_rule": "complete block retained only when every paired truth position occurs in frozen v1 eligible secondary-pair universe; otherwise whole block OOS",
        "formal_blocks_total": len(block_rows),
        "compatible_blocks": len(compatible),
        "oos_blocks": len(out_of_scope),
        "formal_endpoints_total": len(endpoints),
        "compatible_sequences": len(qcache),
        "oos": out_of_scope,
        "sequence_audit": sequence_audit,
    }
    result = {
        "comparator_id": compatibility["comparator_id"],
        "engine": {
            "name": "ViennaRNA Python bindings",
            "version": __import__("RNA").__version__,
            "temperature_c": 37.0,
            "dangles": 2,
            "partition_function": True,
            "constraints": "FORBIDDEN",
            "chemical_reactivity": "FORBIDDEN",
        },
        "population": {
            "formal_endpoint_count": len(endpoints),
            "formal_block_count": len(block_rows),
            "compatible_block_count": len(compatible),
            "oos_block_count": len(out_of_scope),
            "scored_class_contrasts": len(records),
            "scored_parent_count": len({record["member"]["parent_id"] for record in records}),
            "scored_outer_cluster_count": len(cluster_values),
            "compatible_population_identity_sha256": population_identity,
        },
        "result": {
            "F": score,
            "F_95_CI": confidence_interval,
            "bootstrap_unit": "dependency_cluster_F",
            "bootstrap_replicates": 20000,
            "bootstrap_seed": 20261004,
        },
        "provenance": {
            "endpoint_ledger_sha256": sha_file(args.endpoint_ledger),
            "pair_map_sha256": sha_file(args.pair_maps),
            "comparator_code_sha256": sha_file(COMPARATOR_DIR / "vienna_partition_pairing_propensity.py"),
            "frozen_evaluator_code_sha256": sha_file(CODE_DIR / "frozen_aggregation_core.py"),
            "aggregation_trace_sha256": sha_obj(trace),
        },
    }
    (args.output_dir / "RNA05_VIENNA_COMPATIBILITY_AUDIT.json").write_text(
        json.dumps(compatibility, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "RNA05_VIENNA_SCORE_LEDGER.json").write_text(
        json.dumps({"records": records}, separators=(",", ":")) + "\n", encoding="utf-8"
    )
    result["formal_interface"] = "RNA05_F_ONLY_v0.1.3"
    result["score_transform_version"] = "direct_F_credit_v1"
    (args.output_dir / "FORMAL_RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result["result"], sort_keys=True))


if __name__ == "__main__":
    main()
