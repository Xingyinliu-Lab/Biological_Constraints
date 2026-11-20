#!/usr/bin/env python3
"""Single public PS09 entry: model representations to formal scores."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from representations_to_leaf import (
    crossfit,
    endpoint_id,
    load_representations,
    read_folds,
    read_membership,
)
from score_ps09_formal import score, sha256
from ps09_refset import validate_refset


SCORING_DIR = Path(__file__).resolve().parent
REFSET_DIR = SCORING_DIR.parent / "refset"
NUMERIC_BINDING = SCORING_DIR / "PS09_NUMERIC_INPUT_BINDING.json"


def file_record(path: Path) -> dict[str, object]:
    return {
        "path": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the frozen PS09 cross-fit and formal aggregation."
    )
    parser.add_argument("--representations", required=True, type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--refset", type=Path, default=REFSET_DIR)
    args = parser.parse_args()

    paths = {
        "binding": NUMERIC_BINDING,
        "membership": args.refset / "FORMAL_BLOCK_MEMBERSHIP.tsv.gz",
        "blocks": args.refset / "MATCHED_BLOCKS.tsv",
        "clusters": args.refset / "MMSEQS2_CLUSTERS.tsv",
        "fasta": args.refset / "SELECTED_PROTEINS.faa",
        "folds": args.refset / "PS09_FROZEN_CLUSTER_FOLD_ASSIGNMENT.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing PS09 RefSet assets: {missing}")

    closure = validate_refset(
        membership=paths["membership"],
        blocks=paths["blocks"],
        clusters=paths["clusters"],
        fasta=paths["fasta"],
        folds=paths["folds"],
    )
    membership_rows = closure["membership_rows"]
    folds = closure["folds"]
    expected_ids = {endpoint_id(row) for row in membership_rows}
    vectors = load_representations(args.representations, expected_ids)
    heldout_scores = crossfit(membership_rows, folds, vectors)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    leaf_path = args.output_dir / "HELDOUT_ENDPOINT_SCORES.tsv"
    with leaf_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["endpoint_id", "score", "fold"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(
            {
                "endpoint_id": endpoint_id(row),
                "score": repr(float(heldout_scores[index])),
                "fold": folds[row["homology_cluster_id"]],
            }
            for index, row in enumerate(membership_rows)
        )

    summary_path = args.output_dir / "SCORE_SUMMARY.json"
    individual_path = args.output_dir / "INDIVIDUAL_SCORES.tsv"
    pair_path = args.output_dir / "PAIR_CREDIT_LEDGER.tsv.gz"
    block_path = args.output_dir / "BLOCK_SCORES.tsv"
    cluster_path = args.output_dir / "CLUSTER_SCORES.tsv"
    result = score(
        binding=paths["binding"],
        membership=paths["membership"],
        blocks=paths["blocks"],
        clusters=paths["clusters"],
        fasta=paths["fasta"],
        folds=paths["folds"],
        leaf_scores=leaf_path,
        output=summary_path,
        individual_scores=individual_path,
        pair_credits=pair_path,
        block_scores=block_path,
        cluster_scores=cluster_path,
    )
    result["model_id"] = args.model_id
    summary_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    outputs = [leaf_path, individual_path, pair_path, block_path, cluster_path, summary_path]
    bootstrap_path = Path(str(summary_path) + ".bootstrap.float64le.bin")
    outputs.append(bootstrap_path)
    manifest = {
        "constraint_id": "PS09",
        "model_id": args.model_id,
        "formal_interface": "F_ONLY",
        "score_transform_version": "PS09_v0.1.2_F_ONLY",
        "public_entry": Path(__file__).name,
        "representations": {
            "bytes": args.representations.stat().st_size,
            "sha256": sha256(args.representations),
        },
        "outputs": [file_record(path) for path in outputs],
        "refset_closure": closure["counts"],
    }
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
