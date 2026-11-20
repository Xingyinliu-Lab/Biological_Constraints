#!/usr/bin/env python3
"""Extract sequence-free residue embeddings from the ProteinMPNN CA encoder."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ca_from_pdb(path: Path) -> tuple[np.ndarray, list[str]]:
    coordinates: list[list[float]] = []
    residue_keys: list[str] = []
    seen: set[tuple[str, str, str]] = set()
    for line in path.open(errors="replace"):
        if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[16] not in (" ", "A"):
            continue
        key = (line[21], line[22:26].strip(), line[26])
        if key in seen:
            continue
        seen.add(key)
        coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        residue_keys.append(f"{key[0]}:{key[1]}:{key[2].strip()}")
    if not coordinates:
        raise ValueError(f"no CA coordinates: {path}")
    return np.asarray(coordinates, dtype=np.float32), residue_keys


def encode(model, gather_nodes, coordinates: np.ndarray, device: torch.device) -> np.ndarray:
    x = torch.from_numpy(coordinates).unsqueeze(0).to(device)
    length = x.shape[1]
    mask = torch.ones((1, length), dtype=torch.float32, device=device)
    residue_idx = torch.arange(length, device=device, dtype=torch.long).unsqueeze(0)
    chain_encoding = torch.ones((1, length), dtype=torch.long, device=device)
    edges, edge_index = model.features(x, mask, residue_idx, chain_encoding)
    nodes = torch.zeros((edges.shape[0], edges.shape[1], edges.shape[-1]), device=device)
    edges = model.W_e(edges)
    mask_attend = gather_nodes(mask.unsqueeze(-1), edge_index).squeeze(-1)
    mask_attend = mask.unsqueeze(-1) * mask_attend
    for layer in model.encoder_layers:
        nodes, edges = layer(nodes, edges, edge_index, mask, mask_attend)
    return nodes[0].float().cpu().numpy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--proteinmpnn-repo", required=True, type=Path)
    parser.add_argument("--weights", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    args = parser.parse_args()

    sys.path.insert(0, str(args.proteinmpnn_repo))
    from protein_mpnn_utils import ProteinMPNN, gather_nodes

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(args.weights, map_location=device, weights_only=False)
    model = ProteinMPNN(
        ca_only=True,
        num_letters=21,
        node_features=128,
        edge_features=128,
        hidden_dim=128,
        num_encoder_layers=3,
        num_decoder_layers=3,
        augment_eps=0.0,
        k_neighbors=checkpoint["num_edges"],
    ).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    torch.manual_seed(20260914)
    np.random.seed(20260914)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdbs = sorted(args.input_dir.glob("*.pdb"))
    records = []
    started = time.time()
    with torch.inference_mode():
        for number, pdb in enumerate(pdbs, 1):
            output = args.output_dir / f"{pdb.stem}.npz"
            coordinates, residue_keys = ca_from_pdb(pdb)
            if output.is_file():
                with np.load(output, allow_pickle=False) as saved:
                    representation = saved["residue_embedding"]
                if representation.shape != (len(coordinates), 128):
                    raise ValueError(f"stale or malformed checkpoint: {output}")
            else:
                representation = encode(model, gather_nodes, coordinates, device)
                np.savez_compressed(
                    output,
                    residue_embedding=representation,
                    ca_coords=coordinates,
                    residue_keys=np.asarray(residue_keys),
                )
            records.append({
                "object_id": pdb.stem,
                "residue_count": len(coordinates),
                "embedding_dim": int(representation.shape[1]),
                "npz_sha256": sha256(output),
            })
            if number % 50 == 0 or number == len(pdbs):
                print(json.dumps({"done": number, "total": len(pdbs)}), flush=True)

    summary = {
        "status": "PASS",
        "model": "ProteinMPNN-CA-v_48_020 encoder",
        "input_modality": "CA_backbone_coordinates_only",
        "decoder_used": False,
        "sequence_features_used": False,
        "pooling": "none; residue embeddings retained",
        "object_count": len(records),
        "embedding_dim": 128,
        "weights_sha256": sha256(args.weights),
        "source_commit": "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57",
        "elapsed_seconds": time.time() - started,
        "records": records,
    }
    (args.output_dir / "encoding_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "records"}, sort_keys=True))


if __name__ == "__main__":
    main()
