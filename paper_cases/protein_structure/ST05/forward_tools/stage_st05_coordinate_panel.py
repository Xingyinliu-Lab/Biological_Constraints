#!/usr/bin/env python3
"""Flatten the released ST05 base/variant coordinate panel for model encoders."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refset-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    with (args.refset_dir / "objects.tsv").open(encoding="utf-8", newline="") as handle:
        objects = list(csv.DictReader(handle, delimiter="\t"))
    with (args.refset_dir / "relations.tsv").open(encoding="utf-8", newline="") as handle:
        relations = list(csv.DictReader(handle, delimiter="\t"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for row in objects:
        source = args.refset_dir / row["coordinate_locator"]
        shutil.copyfile(source, args.output_dir / f"{row['base_object_id']}__BASE.pdb")
    seen: set[str] = set()
    for row in relations:
        variant_id = row["variant_object_id"]
        if variant_id in seen:
            raise ValueError(f"duplicate variant_object_id: {variant_id}")
        seen.add(variant_id)
        source = args.refset_dir / row["variant_coordinate_locator"]
        shutil.copyfile(source, args.output_dir / f"{variant_id}.pdb")
    print(f"PASS: staged {len(objects)} bases and {len(seen)} variants")


if __name__ == "__main__":
    main()
