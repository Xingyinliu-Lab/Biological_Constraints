from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-key", required=True)
    parser.add_argument("--forward-dir", type=Path, required=True)
    parser.add_argument("--sequence-objects", type=Path, required=True)
    parser.add_argument("--forward-object-table", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("output already exists")
    raw_manifest = json.loads((args.forward_dir / "OUTPUT_MANIFEST.json").read_text(encoding="utf-8"))
    raw_interface = json.loads((args.forward_dir / "MODEL_INTERFACE.json").read_text(encoding="utf-8"))
    if raw_manifest.get("status") != "PASS" or not (args.forward_dir / "RUN_COMPLETE").is_file():
        raise ValueError("forward is not complete")
    if raw_interface.get("model_key") != args.model_key:
        raise ValueError("forward model identity mismatch")
    source = raw_interface.get("source_sha256", {})
    if source.get("DNA02_02_FORWARD_OBJECT_TABLE.tsv") != sha256(args.forward_object_table) or source.get("SOURCE_FASTA") != sha256(args.fasta):
        raise ValueError("forward source identity mismatch")
    files = raw_manifest.get("files", {})
    for name in ["EMBEDDINGS.npz", "REPEAT_EMBEDDINGS.npz", "MODEL_INTERFACE.json", "FORWARD_VALIDATION.json"]:
        path = args.forward_dir / name
        if name not in files or not path.is_file() or files[name]["sha256"] != sha256(path) or int(files[name]["bytes"]) != path.stat().st_size:
            raise ValueError(f"forward manifest mismatch for {name}")
    result = {
        "artifact": "DNA02_02_FORWARD_MANIFEST",
        "version": "0.1.1",
        "status": "PASS",
        "model_key": args.model_key,
        "source_sha256": {
            "sequence_objects": sha256(args.sequence_objects),
            "forward_object_table": sha256(args.forward_object_table),
            "fasta": sha256(args.fasta),
        },
        "raw_forward_interface_sha256": sha256(args.forward_dir / "MODEL_INTERFACE.json"),
        "raw_output_manifest_sha256": sha256(args.forward_dir / "OUTPUT_MANIFEST.json"),
        "files": {name: {"sha256": sha256(args.forward_dir / name), "bytes": (args.forward_dir / name).stat().st_size} for name in files},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=args.output.name + ".staging.", dir=args.output.parent)
    os.close(fd)
    staging = Path(name)
    staging.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(staging, args.output)
    print(json.dumps({"status": "PASS", "model_key": args.model_key, "output": str(args.output)}))


if __name__ == "__main__":
    main()
