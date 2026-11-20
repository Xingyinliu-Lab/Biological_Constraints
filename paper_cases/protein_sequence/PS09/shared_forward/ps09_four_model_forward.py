#!/usr/bin/env python3
"""Configurable PS09 residue Forward for the four registered non-ESM3 models."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from safetensors.numpy import save_file


MODEL_KEYS = {"esm2-650m", "prott5", "esmc-600m", "carp-640m"}


def sha256_path(path: Path) -> str:
    if path.is_file():
        return hashlib.sha256(path.read_bytes()).hexdigest()
    if not path.is_dir():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    for member in sorted(item for item in path.rglob("*") if item.is_file()):
        digest.update(member.relative_to(path).as_posix().encode())
        digest.update(b"\0")
        digest.update(hashlib.sha256(member.read_bytes()).digest())
    return digest.hexdigest()


def window_starts(length: int, width: int = 1022, stride: int = 766) -> list[int]:
    if length <= width:
        return [0]
    starts = [0]
    while starts[-1] + width < length:
        following = min(starts[-1] + stride, length - width)
        if following == starts[-1]:
            break
        starts.append(following)
    return starts


def batches(sequences: dict[str, str], residue_budget: int):
    batch: list[str] = []
    total = 0
    for identity, sequence in sequences.items():
        if batch and total + len(sequence) > residue_budget:
            yield batch
            batch = []
            total = 0
        batch.append(identity)
        total += len(sequence)
    if batch:
        yield batch


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate the registered PS09 residue representations without author-machine paths."
    )
    parser.add_argument("--model", required=True, choices=sorted(MODEL_KEYS))
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--refset", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument(
        "--esm2-import-path",
        type=Path,
        help="Optional directory containing the registered fair-esm Python package.",
    )
    args = parser.parse_args()

    git_root = Path(__file__).resolve().parents[4]
    scoring_dir = git_root / "constraints" / "protein_sequence" / "PS09" / "scoring"
    sys.path.insert(0, str(scoring_dir))
    from ps09_refset import endpoint_id, read_fasta, validate_refset

    paths = {
        "membership": args.refset / "FORMAL_BLOCK_MEMBERSHIP.tsv.gz",
        "blocks": args.refset / "MATCHED_BLOCKS.tsv",
        "clusters": args.refset / "MMSEQS2_CLUSTERS.tsv",
        "fasta": args.refset / "SELECTED_PROTEINS.faa",
        "folds": args.refset / "PS09_FROZEN_CLUSTER_FOLD_ASSIGNMENT.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing PS09 RefSet assets: {missing}")
    closure = validate_refset(**paths)
    rows = closure["membership_rows"]
    sequences = read_fasta(paths["fasta"])
    positions: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for row in rows:
        positions[row["uniprot_accession"]].append(
            (int(row["uniprot_position_1based"]) - 1, endpoint_id(row))
        )

    import torch

    device = torch.device(args.device)
    if device.type == "cuda":
        torch.cuda.set_device(device)
        torch.cuda.reset_peak_memory_stats(device)
    torch.manual_seed(20261008)

    model_metadata: dict[str, object]
    if args.model == "esm2-650m":
        if args.esm2_import_path is not None:
            sys.path.insert(0, str(args.esm2_import_path))
        import esm

        torch.serialization.add_safe_globals([argparse.Namespace])
        network, alphabet = esm.pretrained.load_model_and_alphabet_local(
            str(args.checkpoint)
        )
        network = network.eval().half().to(device)
        converter = alphabet.get_batch_converter()

        def forward(identities: list[str]) -> dict[str, torch.Tensor]:
            output = {}
            for identity in identities:
                sequence = sequences[identity]
                accumulator = None
                counts = None
                for start in window_starts(len(sequence)):
                    part = sequence[start : start + 1022]
                    _, _, tokens = converter([(identity, part)])
                    with torch.inference_mode():
                        hidden = network(
                            tokens.to(device), repr_layers=[33], return_contacts=False
                        )["representations"][33][0, 1 : len(part) + 1].float().cpu()
                    if accumulator is None:
                        accumulator = torch.zeros((len(sequence), hidden.shape[1]))
                        counts = torch.zeros(len(sequence))
                    accumulator[start : start + len(part)] += hidden
                    counts[start : start + len(part)] += 1
                if accumulator is None or counts is None or (counts == 0).any():
                    raise RuntimeError(f"uncovered residue: {identity}")
                output[identity] = accumulator / counts[:, None]
            return output

        iterator = batches(sequences, 512)
        model_metadata = {
            "model": "esm2_t33_650M_UR50D",
            "layer": 33,
            "window_residues": 1022,
            "stride_residues": 766,
            "overlap_rule": "coordinate-anchored mean",
        }
    elif args.model == "prott5":
        from transformers import T5EncoderModel, T5Tokenizer

        tokenizer = T5Tokenizer.from_pretrained(
            str(args.checkpoint), do_lower_case=False, legacy=True, local_files_only=True
        )
        network = T5EncoderModel.from_pretrained(
            str(args.checkpoint), torch_dtype=torch.float16, local_files_only=True
        ).to(device).eval()

        def forward(identities: list[str]) -> dict[str, torch.Tensor]:
            output = {}
            for identity in identities:
                sequence = sequences[identity]
                accumulator = None
                counts = None
                for start in window_starts(len(sequence)):
                    part = sequence[start : start + 1022]
                    encoded = tokenizer(
                        [" ".join(re.sub(r"[UZOB]", "X", part))],
                        add_special_tokens=True,
                        padding=True,
                        return_tensors="pt",
                    )
                    with torch.inference_mode(), torch.autocast(
                        device_type=device.type, dtype=torch.float16, enabled=device.type == "cuda"
                    ):
                        hidden = network(
                            input_ids=encoded["input_ids"].to(device),
                            attention_mask=encoded["attention_mask"].to(device),
                        ).last_hidden_state[0, : len(part)].float().cpu()
                    if accumulator is None:
                        accumulator = torch.zeros((len(sequence), hidden.shape[1]))
                        counts = torch.zeros(len(sequence))
                    accumulator[start : start + len(part)] += hidden
                    counts[start : start + len(part)] += 1
                if accumulator is None or counts is None or (counts == 0).any():
                    raise RuntimeError(f"uncovered residue: {identity}")
                output[identity] = accumulator / counts[:, None]
            return output

        iterator = batches(sequences, 512)
        model_metadata = {
            "model": "Rostlab/prot_t5_xl_half_uniref50-enc",
            "window_residues": 1022,
            "stride_residues": 766,
            "overlap_rule": "coordinate-anchored mean",
        }
    elif args.model == "esmc-600m":
        from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer

        tokenizer = EsmcTokenizer()
        network = EsmcForMaskedLM.from_pretrained(
            str(args.checkpoint), local_files_only=True
        ).to(device).eval()

        def forward(identities: list[str]) -> dict[str, torch.Tensor]:
            encoded = {
                key: value.to(device)
                for key, value in tokenizer(
                    [sequences[identity] for identity in identities],
                    return_tensors="pt",
                    padding=True,
                ).items()
            }
            with torch.inference_mode(), torch.autocast(
                device_type=device.type, dtype=torch.float16, enabled=device.type == "cuda"
            ):
                hidden = network(**encoded).last_hidden_state
            return {
                identity: hidden[index, 1 : len(sequences[identity]) + 1].float().cpu()
                for index, identity in enumerate(identities)
            }

        iterator = batches(sequences, 4096)
        model_metadata = {"model": "biohub/ESMC-600M"}
    else:
        from sequence_models.pretrained import load_model_and_alphabet

        checkpoint = args.checkpoint.resolve()
        if checkpoint.name != "carp_640M.pt":
            raise ValueError("CARP checkpoint must be the registered carp_640M.pt file")
        os.environ["TORCH_HOME"] = str(checkpoint.parents[2])
        network, collater = load_model_and_alphabet("carp_640M")
        network = network.to(device).eval()

        def forward(identities: list[str]) -> dict[str, torch.Tensor]:
            encoded = collater([(sequences[identity],) for identity in identities])[0].to(device)
            with torch.inference_mode():
                hidden = network(encoded, repr_layers=[-1], logits=False)["representations"][
                    len(network.model.embedder.layers)
                ]
            return {
                identity: hidden[index, : len(sequences[identity])].float().cpu()
                for index, identity in enumerate(identities)
            }

        iterator = batches(sequences, 2048)
        model_metadata = {"model": "CARP-640M", "layer": "final_56"}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    destination = args.output_dir / "endpoint_embeddings.safetensors"
    if destination.exists():
        raise FileExistsError(destination)
    started = time.time()
    tensors: dict[str, np.ndarray] = {}
    index_rows: list[dict[str, object]] = []
    for identities in iterator:
        hidden_by_protein = forward(identities)
        for accession, hidden in hidden_by_protein.items():
            if hidden.shape[0] != len(sequences[accession]):
                raise RuntimeError(f"residue length mismatch: {accession}")
            for position, identity in positions[accession]:
                vector = hidden[position].numpy().astype(np.float32)
                if vector.ndim != 1 or vector.size == 0 or not np.isfinite(vector).all():
                    raise RuntimeError(f"invalid endpoint vector: {identity}")
                tensors[identity] = vector
                index_rows.append(
                    {
                        "endpoint_id": identity,
                        "uniprot_accession": accession,
                        "position_0based": position,
                        "embedding_key": identity,
                        "embedding_dim": vector.size,
                    }
                )
    expected = {endpoint_id(row) for row in rows}
    if set(tensors) != expected:
        raise RuntimeError("endpoint closure failure")
    save_file(tensors, str(destination))
    index_path = args.output_dir / "embedding_index.tsv"
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(index_rows)
    report = {
        "status": "PASS",
        "model_id": args.model,
        "formal_endpoint_count": len(index_rows),
        "protein_count": len(sequences),
        "embedding_dimension": len(next(iter(tensors.values()))),
        "checkpoint_sha256": sha256_path(args.checkpoint),
        "embedding_sha256": sha256_path(destination),
        "embedding_index_sha256": sha256_path(index_path),
        "formal_fasta_sha256": sha256_path(paths["fasta"]),
        "seconds": time.time() - started,
        "peak_gpu_mib": (
            torch.cuda.max_memory_allocated(device) / 2**20 if device.type == "cuda" else None
        ),
        "mapping": {
            "coordinate_system": "one-based protein residue from FORMAL_BLOCK_MEMBERSHIP.tsv.gz",
            "special_tokens": "excluded",
            "endpoint_vector": "hidden state at the mapped biological residue",
        },
        **model_metadata,
    }
    (args.output_dir / "FORWARD_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
