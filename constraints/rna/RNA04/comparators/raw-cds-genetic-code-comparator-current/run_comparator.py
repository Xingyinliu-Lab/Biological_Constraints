from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from rna04_raw_cds_genetic_code_comparator import COMPARATOR_ID, score_variant


HERE = Path(__file__).resolve().parent
REFSET = HERE.parents[1] / "refset"
sys.path.insert(0, str(HERE.parents[1] / "scoring"))
from rna04_public_core import aggregate, site_scores, validate_refset


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the complete RNA04 raw-CDS reference")
    parser.add_argument("--ledger", type=Path, default=REFSET / "COMPLETE_SITE_ALT_LEDGER.tsv")
    parser.add_argument("--fasta", type=Path, default=REFSET / "RNA04_unique_endpoints.fasta")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.ledger.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    refset_identity = validate_refset(rows, args.fasta)
    sequences = refset_identity["sequences"]

    endpoint_records: list[dict] = []
    endpoint_scores: list[float] = []
    for row in rows:
        wt = sequences[row["wt_sequence_sha256"]]
        alt = sequences[row["alt_sequence_sha256"]]
        if alt != row["alt_sequence"]:
            raise ValueError(f"ALT sequence mismatch for {row['site_id']}")
        endpoint_id = hashlib.sha256((row["site_id"] + ":" + row["alt_sequence_sha256"]).encode()).hexdigest()
        measured = score_variant(wt, alt)
        endpoint_records.append({"endpoint_id": endpoint_id, "site_id": row["site_id"], "mutation_class": row["mutation_class"], **measured})
        endpoint_scores.append(float(measured["score"]))

    site_values, individual_rows, pair_records = site_scores(rows, np.asarray(endpoint_scores))
    formal_F, clusters, parents = aggregate(rows, site_values)

    args.output.mkdir(parents=True, exist_ok=True)
    endpoint_path = args.output / "ENDPOINT_SCORES.json"
    endpoint_path.write_text(json.dumps({"comparator_id": COMPARATOR_ID, "records": endpoint_records}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_tsv(args.output / "INDIVIDUAL_SCORES.tsv", individual_rows)
    write_tsv(args.output / "PAIR_CREDIT_LEDGER.tsv", pair_records)
    write_tsv(args.output / "PARENT_SCORE_LEDGER.tsv", [{"parent_object_id": key, "parent_F": format(parents[key], ".17g")} for key in sorted(parents)])
    write_tsv(args.output / "CLUSTER_SCORE_LEDGER.tsv", [{"outer_cluster_id": key, "cluster_F": format(clusters[key], ".17g")} for key in sorted(clusters)])
    result = {
        "status": "PASS", "formal_result": True, "assay_id": "C-RNA-04", "comparator_id": COMPARATOR_ID,
        "formal_interface": "F_ONLY", "score_transform_version": "RNA04_v0.2.2_F_ONLY",
        "F": formal_F, "F_95_CI": [formal_F, formal_F],
        "rows": len(rows), "sites": len(site_values), "parents": len(parents), "outer_clusters": len(clusters),
        "ci_method": "degenerate exact cluster-level interval",
        "ledger_sha256": sha256(args.ledger), "fasta_sha256": sha256(args.fasta),
        "population_identity_sha256": refset_identity["population_identity_sha256"],
    }
    (args.output / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    outputs = sorted(path for path in args.output.iterdir() if path.name != "OUTPUT_MANIFEST.json")
    (args.output / "OUTPUT_MANIFEST.json").write_text(json.dumps({
        "files": [{"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs]
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
