"""Run the public frozen PS01 scoring protocol on whole-sequence representations."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from ps01_scoring_core import evaluate


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("invalid or zero-norm representation")
    return float(left @ right / denominator)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refset-dir", type=Path, default=Path(__file__).resolve().parents[1] / "refset")
    parser.add_argument("--representations", required=True, type=Path)
    parser.add_argument("--representation-index", type=Path, help="embedding_index.tsv for executed safetensors Forward output")
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)

    if args.representations.suffix.lower() == ".npz":
        if args.representation_index:
            raise ValueError("--representation-index is only used with safetensors")
        with np.load(args.representations, allow_pickle=False) as archive:
            if set(archive.files) != {"object_ids", "vectors"}:
                raise ValueError("NPZ must contain exactly object_ids and vectors")
            object_ids = [str(value) for value in archive["object_ids"]]
            matrix = np.asarray(archive["vectors"], dtype=np.float64)
    elif args.representations.suffix.lower() == ".safetensors":
        if not args.representation_index:
            raise ValueError("safetensors input requires --representation-index")
        try:
            from safetensors import safe_open
        except ImportError as error:
            raise RuntimeError("install safetensors to read executed paper Forward output") from error
        with args.representation_index.open(encoding="utf-8", newline="") as handle:
            index_rows = list(csv.DictReader(handle, delimiter="\t"))
        if not index_rows or len({row["occurrence_id"] for row in index_rows}) != len(index_rows):
            raise ValueError("invalid or duplicate occurrence IDs in representation index")
        if len({row["embedding_key"] for row in index_rows}) != len(index_rows):
            raise ValueError("duplicate embedding keys in representation index")
        with safe_open(args.representations, framework="numpy") as archive:
            if set(archive.keys()) != {row["embedding_key"] for row in index_rows}:
                raise ValueError("safetensors keys do not match representation index")
            object_ids = [row["occurrence_id"] for row in index_rows]
            matrix = np.stack([np.asarray(archive.get_tensor(row["embedding_key"]), dtype=np.float64) for row in index_rows])
    else:
        raise ValueError("representations must be NPZ or safetensors")
    if (
        matrix.ndim != 2
        or matrix.shape[1] <= 0
        or matrix.shape[0] != len(object_ids)
        or len(object_ids) != len(set(object_ids))
        or not np.isfinite(matrix).all()
    ):
        raise ValueError("invalid representation matrix or object IDs")
    vectors = {object_id: matrix[index] for index, object_id in enumerate(object_ids)}

    relation_path = args.refset_dir / "FORMAL_CONTRAST_LEDGER.tsv.gz"
    with gzip.open(relation_path, "rt", newline="") as handle:
        relations = list(csv.DictReader(handle, delimiter="\t"))
    needed = {value for row in relations for value in (row["a"], row["p"], row["n"])}
    if set(vectors) != needed:
        raise ValueError("representation object universe mismatch")

    ledger_path = args.output_dir / "INDIVIDUAL_SCORES.tsv.gz"
    fields = [
        "clan_id", "anchor_family", "anchor", "positive", "negative_family",
        "negative", "positive_cosine", "negative_cosine", "F_credit",
    ]
    with gzip.open(ledger_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for relation in relations:
            positive = cosine(vectors[relation["a"]], vectors[relation["p"]])
            negative = cosine(vectors[relation["a"]], vectors[relation["n"]])
            writer.writerow(
                {
                    "clan_id": relation["clan_id"],
                    "anchor_family": relation["anchor_family"],
                    "anchor": relation["a"],
                    "positive": relation["p"],
                    "negative_family": relation["negative_family"],
                    "negative": relation["n"],
                    "positive_cosine": repr(positive),
                    "negative_cosine": repr(negative),
                    "F_credit": repr(1.0 if positive > negative else 0.0 if positive == negative else -1.0),
                }
            )

    result_path = args.output_dir / "EVALUATION_RESULT.json"
    result = evaluate(ledger_path, return_bootstrap=True)
    F_bootstrap = np.asarray(result.pop("F_bootstrap"), dtype=np.float64)
    result.update(
        {
            "status": "PASS",
            "evaluation_class": "FORMAL_PROTOCOL_EVALUATION",
            "formal_scoring_protocol": True,
            "registered_paper_result": False,
            "constraint_id": "PS01",
            "model_id": args.model_id,
            "coverage": {
                "formal_contrasts": len(relations),
                "formal_objects": len(needed),
                "scored_contrasts": len(relations),
                "missing_objects": 0,
            },
            "inputs": {
                "refset_release": "PS01_ASSAY_REFSET_RELEASE_v2.0.0_BROAD_16K_20261006",
                "relation_ledger_sha256": hashlib.sha256(relation_path.read_bytes()).hexdigest(),
                "representations_sha256": hashlib.sha256(args.representations.read_bytes()).hexdigest(),
                "representation_index_sha256": hashlib.sha256(args.representation_index.read_bytes()).hexdigest() if args.representation_index else None,
            },
        }
    )
    np.save(args.output_dir / "F_BOOTSTRAP.npy", F_bootstrap, allow_pickle=False)
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    files = [ledger_path, args.output_dir / "F_BOOTSTRAP.npy", result_path]
    manifest = {
        "status": "PASS",
        "files": [
            {"path": path.name, "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in files
        ]
    }
    (args.output_dir / "OUTPUT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
