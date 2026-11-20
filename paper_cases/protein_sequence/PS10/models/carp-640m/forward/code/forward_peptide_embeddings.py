#!/usr/bin/env python3
"""Produce mean-pooled peptide embeddings for the frozen PS10 population.

This driver intentionally writes only record identifiers and numeric embeddings.
"""
import argparse
import csv
import hashlib
import json
import os
import pathlib
import re
import sys
import time

import numpy as np


EXPECTED_RECORDS_SHA256 = "c988e89f50d92ecdda8805ba82ae7f9f95483eea5bd342bdcedfd7e41c9a0dff"
EXPECTED_COUNT = 2031
MODELS = ("esm2", "prott5", "esmc", "carp640")


def sha256(path: pathlib.Path) -> str:
    if path.is_dir():
        h = hashlib.sha256()
        for child in sorted((p for p in path.rglob("*") if p.is_file()), key=lambda p: p.relative_to(path).as_posix()):
            h.update(child.relative_to(path).as_posix().encode("utf-8"))
            h.update(b"\0")
            h.update(sha256(child).encode("ascii"))
            h.update(b"\n")
        return h.hexdigest()
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def batches(rows, limit):
    group, used = [], 0
    for row in rows:
        width = len(row["normalized_peptide_identity"])
        if group and used + width > limit:
            yield group
            group, used = [], 0
        group.append(row)
        used += width
    if group:
        yield group


def load_rows(path):
    if sha256(path) != EXPECTED_RECORDS_SHA256:
        raise RuntimeError("frozen records SHA256 mismatch")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row["record_id"] for row in rows]
    if len(rows) != EXPECTED_COUNT or len(set(ids)) != EXPECTED_COUNT:
        raise RuntimeError("frozen record population mismatch")
    if any(not re.fullmatch(r"[A-Z]+", row["normalized_peptide_identity"]) for row in rows):
        raise RuntimeError("invalid normalized peptide identity")
    return rows


def encode(model_name, checkpoint, rows, torch_home=None):
    import torch

    torch.cuda.set_device(0)
    torch.cuda.reset_peak_memory_stats()
    torch.manual_seed(20261008)
    if model_name == "esm2":
        import esm

        torch.serialization.add_safe_globals([argparse.Namespace])
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(checkpoint))
        model = model.eval().half().cuda()
        convert = alphabet.get_batch_converter()

        def forward(chunk):
            _, _, tokens = convert([(r["record_id"], r["normalized_peptide_identity"]) for r in chunk])
            with torch.inference_mode():
                hidden = model(tokens.cuda(), repr_layers=[33], return_contacts=False)["representations"][33]
            return [hidden[i, 1 : len(row["normalized_peptide_identity"]) + 1].float().cpu() for i, row in enumerate(chunk)]

        meta = {"model": "esm2_t33_650M_UR50D", "layer": 33, "token_limit": 4096}
    elif model_name == "prott5":
        from transformers import T5EncoderModel, T5Tokenizer

        tokenizer = T5Tokenizer.from_pretrained(str(checkpoint), do_lower_case=False, legacy=True, local_files_only=True)
        model = T5EncoderModel.from_pretrained(str(checkpoint), torch_dtype=torch.float16, local_files_only=True).cuda().eval()

        def forward(chunk):
            prepared = [" ".join(re.sub(r"[UZOB]", "X", row["normalized_peptide_identity"])) for row in chunk]
            tokenized = tokenizer(prepared, add_special_tokens=True, padding=True, return_tensors="pt")
            with torch.inference_mode(), torch.autocast("cuda", dtype=torch.float16):
                hidden = model(input_ids=tokenized["input_ids"].cuda(), attention_mask=tokenized["attention_mask"].cuda()).last_hidden_state
            return [hidden[i, : len(row["normalized_peptide_identity"])].float().cpu() for i, row in enumerate(chunk)]

        meta = {"model": "Rostlab/prot_t5_xl_half_uniref50-enc", "token_limit": 4096}
    elif model_name == "esmc":
        from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer

        tokenizer = EsmcTokenizer()
        model = EsmcForMaskedLM.from_pretrained(str(checkpoint), local_files_only=True).cuda().eval()

        def forward(chunk):
            tokenized = {key: value.cuda() for key, value in tokenizer([row["normalized_peptide_identity"] for row in chunk], return_tensors="pt", padding=True).items()}
            with torch.inference_mode(), torch.autocast("cuda", dtype=torch.float16):
                hidden = model(**tokenized).last_hidden_state
            return [hidden[i, 1 : len(row["normalized_peptide_identity"]) + 1].float().cpu() for i, row in enumerate(chunk)]

        meta = {"model": "biohub/ESMC-600M", "token_limit": 4096}
    elif model_name == "carp640":
        from sequence_models.pretrained import load_model_and_alphabet

        torch_home = torch_home or checkpoint.parents[2]
        os.environ["TORCH_HOME"] = str(torch_home)
        expected_cache = pathlib.Path(torch.hub.get_dir()) / "checkpoints" / "carp_640M.pt"
        if expected_cache.resolve() != checkpoint.resolve():
            raise RuntimeError("CARP checkpoint must be TORCH_HOME/hub/checkpoints/carp_640M.pt")
        model, collater = load_model_and_alphabet("carp_640M")
        model = model.cuda().eval()

        def forward(chunk):
            tokens = collater([(row["normalized_peptide_identity"],) for row in chunk])[0].cuda()
            with torch.inference_mode():
                hidden = model(tokens, repr_layers=[-1], logits=False)["representations"][len(model.model.embedder.layers)]
            return [hidden[i, : len(row["normalized_peptide_identity"])].float().cpu() for i, row in enumerate(chunk)]

        meta = {"model": "CARP-640M", "layer": "final_56", "token_limit": 2048}
    else:
        raise ValueError(model_name)

    pooled = []
    for chunk in batches(rows, meta["token_limit"]):
        hidden = forward(chunk)
        for row, values in zip(chunk, hidden, strict=True):
            if values.shape[0] != len(row["normalized_peptide_identity"]):
                raise RuntimeError("residue-length mismatch")
            vector = values.mean(dim=0).numpy().astype(np.float32, copy=False)
            if not np.isfinite(vector).all():
                raise RuntimeError("non-finite pooled embedding")
            pooled.append(vector)
    matrix = np.stack(pooled)
    if matrix.shape[0] != EXPECTED_COUNT or not np.isfinite(matrix).all():
        raise RuntimeError("embedding population closure failure")
    return matrix, meta, float(torch.cuda.max_memory_allocated() / 2**20)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, choices=MODELS)
    parser.add_argument("--checkpoint", required=True, type=pathlib.Path)
    parser.add_argument("--torch-home", type=pathlib.Path, help="CARP cache root; defaults from --checkpoint")
    parser.add_argument("--records", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    if not args.checkpoint.exists():
        raise RuntimeError("checkpoint unavailable")
    rows = load_rows(args.records)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "embeddings.npz"
    report_path = args.output_dir / "FORWARD_REPORT.json"
    if output.exists() or report_path.exists():
        raise RuntimeError("refusing to overwrite existing forward output")
    started = time.time()
    matrix, meta, peak_gpu_mib = encode(args.model, args.checkpoint, rows, args.torch_home)
    temp = args.output_dir / "embeddings.tmp.npz"
    np.savez(temp, record_id=np.array([row["record_id"] for row in rows]), embedding=matrix)
    temp.replace(output)
    report = {
        "status": "PASS",
        "model_key": args.model,
        "frozen_record_count": len(rows),
        "embedding_dim": int(matrix.shape[1]),
        "records_sha256": sha256(args.records),
        "checkpoint_path": str(args.checkpoint.resolve()),
        "checkpoint_sha256": sha256(args.checkpoint),
        "embeddings_sha256": sha256(output),
        "seconds": time.time() - started,
        "peak_gpu_mib": peak_gpu_mib,
        "pooling": "unweighted_mean_over_residue_embeddings",
        **meta,
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
