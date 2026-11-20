#!/usr/bin/env python3
"""Parameterized PS01 whole-protein Forward exporter for four paper cases."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from pathlib import Path

import torch
from safetensors.torch import save_file


def sha256(path: Path) -> str:
    path = Path(path)
    digest = hashlib.sha256()
    if path.is_dir():
        for child in sorted(item for item in path.rglob("*") if item.is_file()):
            digest.update(child.relative_to(path).as_posix().encode("utf-8"))
            digest.update(bytes.fromhex(sha256(child)))
        return digest.hexdigest()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: dict[str, str] = {}
    identifier = None
    sequence: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(sequence)
                identifier, sequence = line[1:].split()[0], []
            elif line:
                sequence.append(line)
    if identifier is not None:
        records[identifier] = "".join(sequence)
    if not records or len(records) != len(set(records)):
        raise ValueError("empty or duplicate-ID FASTA")
    return records


def batches(sequences: dict[str, str], count: int, tokens: int):
    batch: list[str] = []
    used = 0
    for endpoint_id in sorted(sequences, key=lambda key: (len(sequences[key]), key)):
        size = len(sequences[endpoint_id]) + 2
        if batch and (len(batch) >= count or used + size > tokens):
            yield batch
            batch, used = [], 0
        batch.append(endpoint_id)
        used += size
    if batch:
        yield batch


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=("esm2", "prott5", "esmc", "carp"), required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"output directory already exists: {args.output_dir}")
    if not args.checkpoint.exists():
        raise FileNotFoundError(args.checkpoint)
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required")

    sequences = read_fasta(args.fasta)
    if len(sequences) != 2316:
        raise ValueError(f"expected 2316 PS01 objects, found {len(sequences)}")
    torch.cuda.set_device(0)
    torch.cuda.reset_peak_memory_stats()
    if args.model == "esm2":
        import esm

        model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(args.checkpoint))
        model = model.eval().half().cuda()
        converter = alphabet.get_batch_converter()
        meta, iterator = {"pooling": "residue_mean", "layer": "33"}, batches(sequences, 16, 4096)

        def forward(ids: list[str]):
            _, _, token_ids = converter([(key, sequences[key]) for key in ids])
            with torch.inference_mode():
                hidden = model(token_ids.cuda(), repr_layers=[33], return_contacts=False)["representations"][33]
            return [hidden[i, 1 : len(sequences[key]) + 1].float().mean(0).cpu().contiguous() for i, key in enumerate(ids)]
    elif args.model == "prott5":
        from transformers import T5EncoderModel, T5Tokenizer

        tokenizer = T5Tokenizer.from_pretrained(str(args.checkpoint), do_lower_case=False, legacy=True, local_files_only=True)
        model = T5EncoderModel.from_pretrained(str(args.checkpoint), torch_dtype=torch.float16, local_files_only=True).cuda().eval()
        meta, iterator = {"pooling": "residue_mean", "UO_mapping": "X"}, batches(sequences, 8, 2048)

        def forward(ids: list[str]):
            inputs = tokenizer([" ".join(re.sub(r"[UZOB]", "X", sequences[key])) for key in ids], add_special_tokens=True, padding=True, return_tensors="pt")
            with torch.inference_mode(), torch.autocast("cuda", dtype=torch.float16):
                hidden = model(input_ids=inputs["input_ids"].cuda(), attention_mask=inputs["attention_mask"].cuda()).last_hidden_state
            return [hidden[i, : len(sequences[key])].float().mean(0).cpu().contiguous() for i, key in enumerate(ids)]
    elif args.model == "esmc":
        from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer

        tokenizer = EsmcTokenizer()
        model = EsmcForMaskedLM.from_pretrained(str(args.checkpoint), local_files_only=True).cuda().eval()
        meta, iterator = {"pooling": "residue_mean"}, batches(sequences, 16, 4096)

        def forward(ids: list[str]):
            inputs = {key: value.cuda() for key, value in tokenizer([sequences[key] for key in ids], return_tensors="pt", padding=True).items()}
            with torch.inference_mode(), torch.autocast("cuda", dtype=torch.float16):
                hidden = model(**inputs).last_hidden_state
            return [hidden[i, 1 : len(sequences[key]) + 1].float().mean(0).cpu().contiguous() for i, key in enumerate(ids)]
    else:
        from sequence_models.pretrained import load_model_and_alphabet

        if args.checkpoint.name != "carp_640M.pt":
            raise ValueError("CARP checkpoint must be the registered carp_640M.pt file")
        os.environ["TORCH_HOME"] = str(args.checkpoint.parents[2])
        model, collater = load_model_and_alphabet("carp_640M")
        model = model.cuda().eval()
        meta, iterator = {"pooling": "residue_mean", "layer": "final_56"}, batches(sequences, 4, 2048)

        def forward(ids: list[str]):
            token_ids = collater([(sequences[key],) for key in ids])[0].cuda()
            with torch.inference_mode():
                hidden = model(token_ids, repr_layers=[-1], logits=False)["representations"][len(model.model.embedder.layers)]
            return [hidden[i, : len(sequences[key])].float().mean(0).cpu().contiguous() for i, key in enumerate(ids)]

    vectors, index = {}, []
    for ids in iterator:
        for endpoint_id, vector in zip(ids, forward(ids)):
            if not torch.isfinite(vector).all():
                raise RuntimeError(f"non-finite vector: {endpoint_id}")
            key = f"embedding_{len(vectors) + 1:06d}"
            vectors[key] = vector
            index.append((endpoint_id, key, len(sequences[endpoint_id]), vector.numel(), hashlib.sha256(vector.numpy().tobytes()).hexdigest()))
    args.output_dir.mkdir(parents=True)
    embedding_path = args.output_dir / "embeddings.safetensors"
    index_path = args.output_dir / "embedding_index.tsv"
    save_file(vectors, str(embedding_path), metadata=meta)
    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("occurrence_id", "embedding_key", "sequence_length", "embedding_dim", "embedding_sha256"))
        writer.writerows(index)
    report = {"status": "PASS", "constraint": "PS01", "interface_id": "PS01_PUBLIC_WHOLE_PROTEIN_FORWARD_v1.0.0", "model": args.model, "endpoint_count": len(index), "pooling": "residue_mean", "input_fasta_sha256": sha256(args.fasta), "checkpoint_sha256": sha256(args.checkpoint), "embedding_sha256": sha256(embedding_path), "embedding_index_sha256": sha256(index_path)}
    (args.output_dir / "FORWARD_REPORT.json").write_text(__import__("json").dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

