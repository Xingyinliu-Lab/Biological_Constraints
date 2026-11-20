#!/usr/bin/env python3
"""Convert model representations into the frozen PS08 leaf-score ledger."""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from ps08_contract import ID_FIELDS, load_population, sha256

SCORE_FIELDS = (
    "left_near_score",
    "left_far_score",
    "right_near_score",
    "right_far_score",
)
def load_embeddings(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {"object_ids", "embeddings"}:
            raise ValueError("NPZ must contain exactly object_ids and embeddings")
        object_ids = np.asarray(data["object_ids"])
        raw_embeddings = np.asarray(data["embeddings"])
        if not np.issubdtype(raw_embeddings.dtype, np.number) or np.iscomplexobj(raw_embeddings):
            raise ValueError("embeddings must be real numeric vectors")
        embeddings = np.asarray(raw_embeddings, dtype=np.float64)
    if object_ids.ndim != 1 or embeddings.ndim != 2:
        raise ValueError("object_ids must be 1-D and embeddings must be 2-D")
    if embeddings.shape[1] <= 0:
        raise ValueError("embeddings must have a positive feature dimension")
    if len(object_ids) != len(embeddings) or len(object_ids) == 0:
        raise ValueError("object_ids and embeddings must have the same non-zero row count")
    ids = [str(value) for value in object_ids.tolist()]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate object_id in representation NPZ")
    if not np.isfinite(embeddings).all():
        raise ValueError("representations contain NaN or Inf")
    norms = np.linalg.norm(embeddings, axis=1)
    if not np.isfinite(norms).all():
        raise ValueError("representation norms are non-finite")
    normalized = np.zeros_like(embeddings, dtype=np.float64)
    nonzero = norms > 0
    normalized[nonzero] = embeddings[nonzero] / norms[nonzero, None]
    return dict(zip(ids, normalized, strict=True))


def cosine(vectors: dict[str, np.ndarray], left: str, right: str) -> float:
    try:
        value = float(np.dot(vectors[left], vectors[right]))
    except KeyError as exc:
        raise ValueError(f"missing representation for object_id {exc.args[0]}") from exc
    if not math.isfinite(value):
        raise ValueError("non-finite cosine score")
    return value


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create FORMAL_LEAF_SCORES.tsv from PS08 object representations."
    )
    parser.add_argument("--population", type=Path, required=True)
    parser.add_argument("--objects-fasta", type=Path, required=True)
    parser.add_argument("--representations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    members, population_metadata = load_population(args.population, args.objects_fasta)

    vectors = load_embeddings(args.representations)
    expected_ids = {
        row[field]
        for row in members
        for field in ("x0_object_id", "x1_object_id", "x2_object_id")
    }
    if set(vectors) != expected_ids:
        missing = len(expected_ids - set(vectors))
        extra = len(set(vectors) - expected_ids)
        raise ValueError(f"representation object universe mismatch: missing={missing}, extra={extra}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*ID_FIELDS, *SCORE_FIELDS], delimiter="\t")
        writer.writeheader()
        for row in members:
            x0 = row["x0_object_id"]
            x1 = row["x1_object_id"]
            x2 = row["x2_object_id"]
            writer.writerow(
                {
                    **{field: row[field] for field in ID_FIELDS},
                    "left_near_score": format(cosine(vectors, x0, x1), ".17g"),
                    "left_far_score": format(cosine(vectors, x0, x2), ".17g"),
                    "right_near_score": format(cosine(vectors, x2, x1), ".17g"),
                    "right_far_score": format(cosine(vectors, x2, x0), ".17g"),
                }
            )
    metadata = {
        "kind": "PS08_REPRESENTATIONS_TO_LEAF_METADATA",
        "population": population_metadata,
        "representations_sha256": sha256(args.representations),
        "leaf_ledger_sha256": sha256(args.output),
        "leaf_count": len(members),
        "object_count": len(expected_ids),
    }
    args.output.with_name(args.output.name + ".metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
