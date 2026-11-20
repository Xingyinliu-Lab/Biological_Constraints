#!/usr/bin/env python3
"""Encode an incremental PDB directory with the frozen ESM-IF1 encoder."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch


def coords_from_pdb(path: Path) -> np.ndarray:
    residues, order = {}, []
    for line in path.open(errors="replace"):
        if not line.startswith("ATOM"): continue
        atom = line[12:16].strip()
        if atom not in ("N", "CA", "C") or line[16] not in (" ", "A"): continue
        key = (line[21], line[22:26], line[26])
        if key not in residues: residues[key] = {}; order.append(key)
        if atom not in residues[key]:
            residues[key][atom] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    good = [np.asarray([residues[k][a] for a in ("N", "CA", "C")], dtype=np.float32)
            for k in order if all(a in residues[k] for a in ("N", "CA", "C"))]
    if not good: raise ValueError(f"no complete backbone residues: {path}")
    return np.stack(good)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--esm-repo", required=True, type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.esm_repo))
    import esm
    from esm.inverse_folding.util import get_encoder_output

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdbs = sorted(args.input_dir.glob("*.pdb"))
    torch.manual_seed(20260910); np.random.seed(20260910)
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model.eval().cuda()
    completed, started = [], time.time()
    for number, path in enumerate(pdbs, 1):
        output = args.output_dir / f"{path.stem}.npz"
        if output.is_file():
            with np.load(output) as saved: residue_count = len(saved["residue_embedding"])
        else:
            coordinates = coords_from_pdb(path)
            with torch.inference_mode():
                representation = get_encoder_output(model, alphabet, coordinates).float().cpu().numpy()
            np.savez_compressed(output, residue_embedding=representation, coords=coordinates)
            residue_count = len(representation)
        completed.append({"object_id": path.stem, "residue_count": residue_count,
                          "npz_sha256": hashlib.sha256(output.read_bytes()).hexdigest()})
        if number % 10 == 0: print(json.dumps({"done": number, "total": len(pdbs)}), flush=True)
    summary = {"status": "PASS", "objects": len(completed), "elapsed_seconds": time.time() - started,
               "records": completed}
    (args.output_dir / "incremental_encoding_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "records"}, sort_keys=True))


if __name__ == "__main__":
    main()
