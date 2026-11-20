#!/usr/bin/env python3
"""Stage ST02 coordinate copies for ProToken without altering the released RefSet."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def normalize_insertion_codes(text: str) -> tuple[str, bool]:
    """Sequentially renumber residues within chain only when insertion codes occur."""
    if not any(line.startswith(("ATOM  ", "HETATM")) and line[26:27].strip() for line in text.splitlines()):
        return text, False
    mapping: dict[tuple[str, str, str], int] = {}
    next_number: dict[str, int] = {}
    output: list[str] = []
    for line in text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            output.append(line)
            continue
        chain, residue, insertion = line[21:22], line[22:26], line[26:27]
        key = (chain, residue, insertion)
        if key not in mapping:
            next_number[chain] = next_number.get(chain, 0) + 1
            mapping[key] = next_number[chain]
        output.append(f"{line[:22]}{mapping[key]:4d} {line[27:]}")
    return "".join(output), True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--objects", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    with args.objects.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    normalized: list[str] = []
    for row in rows:
        locator = Path(row["coordinate_locator"])
        source = args.input_dir / locator.name
        if not source.is_file():
            raise FileNotFoundError(source)
        staged, changed = normalize_insertion_codes(source.read_text(encoding="utf-8"))
        (args.output_dir / locator.name).write_text(staged, encoding="utf-8")
        if changed:
            normalized.append(row["record_id"])
    print(json.dumps({
        "status": "PASS", "input_objects": len(rows),
        "normalized_insertion_code_objects": normalized,
        "policy": "copied inputs only; sequential residue renumbering within chain for PDB insertion codes",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
