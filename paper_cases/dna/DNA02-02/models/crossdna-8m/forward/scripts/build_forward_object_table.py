from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-objects", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("output already exists")
    with args.sequence_objects.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"object_id", "sequence_sha256", "sequence_length", "normalized_alphabet"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("unexpected sequence object schema")
    ids: set[str] = set()
    normalized = []
    for row in rows:
        if row["object_id"] in ids or row["normalized_alphabet"] != "ACGT":
            raise ValueError("duplicate object or invalid alphabet")
        ids.add(row["object_id"])
        normalized.append({"object_id": row["object_id"], "length_bp": int(row["sequence_length"]), "sequence_sha256": row["sequence_sha256"]})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fd, staging_name = tempfile.mkstemp(prefix=args.output.name + ".staging.", dir=args.output.parent)
    os.close(fd)
    staging = Path(staging_name)
    with staging.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["object_id", "length_bp", "sequence_sha256"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized)
    os.replace(staging, args.output)


if __name__ == "__main__":
    main()
