from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path


RUNNER_ID = "ST02_TMALIGN_PAIR_RUNNER_v0.1.3_F_ONLY"
TM_PATTERN = re.compile(r"TM-score=\s*([0-9.eE+-]+)\s*\((?:if )?normalized by length of (?:Chain|Structure)_([12])")
HERE = Path(__file__).resolve().parent
CONSTRAINT = HERE.parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_coordinate(objects_path: Path, locator: str) -> Path:
    path = Path(locator)
    if not path.is_absolute():
        path = objects_path.parent / path
    return path.resolve()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmalign", required=True, type=Path, help="Local TM-align executable")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    executable = args.tmalign
    if not executable.is_file():
        raise FileNotFoundError("supply --tmalign with a TM-align executable")
    executable_sha256 = sha256(executable)
    objects_path = CONSTRAINT / "refset" / "objects.tsv"
    relations_path = CONSTRAINT / "refset" / "relations.tsv"
    with objects_path.open(newline="", encoding="utf-8") as handle:
        object_rows = list(csv.DictReader(handle, delimiter="\t"))
    objects = {row["record_id"]: row for row in object_rows}
    if len(objects) != len(object_rows):
        raise ValueError("duplicate object ID")
    with relations_path.open(newline="", encoding="utf-8") as handle:
        relations = list(csv.DictReader(handle, delimiter="\t"))
    needed = set()
    for row in relations:
        needed.add(tuple(sorted((row["anchor_object_id"], row["true_target_object_id"]))))
        needed.add(tuple(sorted((row["anchor_object_id"], row["negative_object_id"]))))

    prefix: list[str] = []
    pair_rows = []
    for object_a, object_b in sorted(needed):
        coordinates = []
        for object_id in (object_a, object_b):
            if object_id not in objects:
                raise ValueError("relation object missing from objects table")
            path = resolve_coordinate(objects_path, objects[object_id]["coordinate_locator"])
            if sha256(path) != objects[object_id]["coordinate_sha256"]:
                raise ValueError("coordinate authority mismatch")
            coordinates.append(path)
        completed = subprocess.run(
            [str(executable), *prefix, str(coordinates[0]), str(coordinates[1])],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False,
        )
        if completed.returncode != 0:
            raise ValueError("TM-align execution failure")
        parsed = {int(chain): float(value) for value, chain in TM_PATTERN.findall(completed.stdout)}
        if set(parsed) != {1, 2} or not all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in parsed.values()):
            raise ValueError("TM-align output parse/nonfinite failure")
        pair_rows.append({
            "object_a": object_a, "object_b": object_b,
            "tm_score_chain_1_normalized": parsed[1],
            "tm_score_chain_2_normalized": parsed[2],
            "symmetric_tm_score": (parsed[1] + parsed[2]) / 2.0,
        })
    result = {
        "artifact_id": "ST02_TMALIGN_PAIR_LEDGER",
        "status": "COMPLETE_PAIR_LEDGER",
        "runner_id": RUNNER_ID,
        "pair_runner_sha256": sha256(Path(__file__).resolve()),
        "tmalign_executable": {
            "path": str(executable.resolve()),
            "sha256": executable_sha256,
        },
        "required_pair_count": len(needed), "pair_count": len(pair_rows),
        "pairs": pair_rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
