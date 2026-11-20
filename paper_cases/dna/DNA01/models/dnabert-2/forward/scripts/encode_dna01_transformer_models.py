#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch

# Use PyTorch math-SDP for the registered standard-attention fallback.
if torch.cuda.is_available():
    torch.backends.cuda.enable_flash_sdp(False)
    torch.backends.cuda.enable_mem_efficient_sdp(False)
    torch.backends.cuda.enable_math_sdp(True)


MODEL_IDS = {
    "nt": "InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
    "dnabert2": "zhihan1996/DNABERT-2-117M",
    "caduceus_ps": "kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
}
BATCH_SIZES = {"nt": 64, "dnabert2": 48, "caduceus_ps": 24}
WINDOW_BP = 1000
DNABERT2_FORMAL_REVISION = "7bce263b15377fc15361f52cfab88f8b586abda0"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path):
    records, name, chunks = {}, None, []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if name is not None:
                if name in records:
                    raise ValueError("duplicate FASTA identifier")
                records[name] = "".join(chunks).upper()
            name, chunks = raw[1:].split()[0], []
        else:
            chunks.append(raw.strip())
    if name is not None:
        if name in records:
            raise ValueError("duplicate FASTA identifier")
        records[name] = "".join(chunks).upper()
    return records


def read_objects(path, fasta):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids, seen = [], set()
    for row in rows:
        oid = row["object_id"]
        if oid in seen or oid not in fasta:
            raise ValueError("invalid object universe")
        sequence = fasta[oid]
        if len(sequence) != int(row["length_bp"]) or hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
            raise ValueError("object sequence identity mismatch")
        if set(sequence) - set("ACGT"):
            raise ValueError("non-ACGT sequence")
        seen.add(oid)
        ids.append(oid)
    return rows, ids


def load_model(name, dnabert2_model_path=None):
    model_id = MODEL_IDS[name]
    if name == "dnabert2":
        if dnabert2_model_path is None:
            raise ValueError("--dnabert2-model-path is required for the formal Direct HF execution")
        from transformers import AutoModel, AutoTokenizer
        model_source = str(dnabert2_model_path)
        tokenizer = AutoTokenizer.from_pretrained(model_source, revision=DNABERT2_FORMAL_REVISION, trust_remote_code=True, local_files_only=True)
        model = AutoModel.from_pretrained(model_source, revision=DNABERT2_FORMAL_REVISION, trust_remote_code=True, local_files_only=True)
        # DNABERT-2 bundles a Triton FlashAttention kernel outside PyTorch's
        # deterministic-algorithm controls.  Its own implementation provides
        # the mathematically equivalent standard-attention fallback.
        model_module = sys.modules[type(model).__module__]
        if hasattr(model_module, "flash_attn_qkvpacked_func"):
            model_module.flash_attn_qkvpacked_func = None
    else:
        from transformers import AutoModelForMaskedLM, AutoTokenizer
        tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True, local_files_only=True)
        model = AutoModelForMaskedLM.from_pretrained(model_id, trust_remote_code=True, local_files_only=True)
    model.eval().cuda()
    if name == "dnabert2":
        model.float()
    else:
        model.half()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model_id, tokenizer, model


def encode_pass(name, tokenizer, model, object_ids, fasta, batch_size):
    torch.manual_seed(20261013)
    torch.cuda.manual_seed_all(20261013)
    tasks, weights = [], defaultdict(int)
    for oid in object_ids:
        sequence = fasta[oid]
        for start in range(0, len(sequence), WINDOW_BP):
            window = sequence[start:start + WINDOW_BP]
            tasks.append((len(window), oid, window))
            weights[oid] += len(window)
    tasks.sort(key=lambda item: (item[0], item[1], item[2]))
    special = set(tokenizer.all_special_ids)
    sums = {}
    for start in range(0, len(tasks), batch_size):
        batch = tasks[start:start + batch_size]
        encoded = tokenizer([item[2] for item in batch], return_tensors="pt", padding=True, truncation=False)
        encoded = {key: value.cuda() for key, value in encoded.items()}
        use_fp16 = name != "dnabert2"
        with torch.inference_mode(), torch.autocast(device_type="cuda", dtype=torch.float16, enabled=use_fp16):
            output = model(**encoded, output_hidden_states=(name != "dnabert2"), return_dict=True)
            hidden = (output.last_hidden_state if hasattr(output, "last_hidden_state") else output[0]) if name == "dnabert2" else output.hidden_states[-1]
        mask = encoded.get("attention_mask", torch.ones(hidden.shape[:2], device=hidden.device, dtype=torch.long)).bool()
        for token_id in special:
            mask &= encoded["input_ids"].ne(token_id)
        pooled = (hidden.float() * mask.unsqueeze(-1)).sum(1) / mask.sum(1, keepdim=True).clamp_min(1)
        for (bp, oid, _), vector in zip(batch, pooled.cpu()):
            sums[oid] = sums.get(oid, 0) + vector * bp
        if start == 0 or start + len(batch) == len(tasks) or (start // batch_size) % 20 == 0:
            print(json.dumps({"model": name, "windows_done": min(start + len(batch), len(tasks)), "windows_total": len(tasks)}), flush=True)
    vectors = torch.stack([(sums[oid] / weights[oid]).float() for oid in object_ids]).numpy()
    if not np.isfinite(vectors).all():
        raise ValueError("non-finite embedding")
    return vectors, len(tasks)


def write_metadata(stage, names):
    files = {name: {"sha256": sha256_file(stage / name), "bytes": (stage / name).stat().st_size} for name in names}
    (stage / "OUTPUT_MANIFEST.json").write_text(json.dumps({"status": "PASS", "files": files}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    names = names + ["OUTPUT_MANIFEST.json"]
    (stage / "SHA256SUMS").write_text("".join(f"{sha256_file(stage / name)}  {name}\n" for name in names), encoding="utf-8")
    (stage / "RUN_COMPLETE").write_text("PASS\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=sorted(MODEL_IDS), required=True)
    parser.add_argument("--object-table", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dnabert2-model-path", type=Path)
    parser.add_argument("--strict-repeat", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("output path must not already exist")
    if args.strict_repeat:
        torch.use_deterministic_algorithms(True)
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
    fasta = read_fasta(args.fasta)
    rows, object_ids = read_objects(args.object_table, fasta)
    model_id, tokenizer, model = load_model(args.model, args.dnabert2_model_path)
    started = time.time()
    primary, window_count = encode_pass(args.model, tokenizer, model, object_ids, fasta, BATCH_SIZES[args.model])
    repeat, repeat_window_count = encode_pass(args.model, tokenizer, model, object_ids, fasta, BATCH_SIZES[args.model])
    max_delta = float(np.max(np.abs(primary - repeat)))
    repeat_exact = bool(np.array_equal(primary, repeat))
    if args.strict_repeat and not repeat_exact:
        raise RuntimeError(f"strict repeat mismatch: max_abs_difference={max_delta}")
    interface = {
        "artifact": "DNA01_MODEL_INTERFACE",
        "status": "EXECUTED_FROZEN_INTERFACE",
        "model_key": args.model,
        "model_id": model_id,
        "model_revision": DNABERT2_FORMAL_REVISION if args.model == "dnabert2" else getattr(model.config, "_commit_hash", None),
        "representation": "last hidden layer token mean within each nonoverlapping 1000-bp window, then bp-weighted window mean",
        "window_bp": WINDOW_BP,
        "special_token_policy": "exclude tokenizer special tokens",
        "precision": "model float32; pooling and released vectors float32" if args.model == "dnabert2" else "model float16; pooling and released vectors float32",
        "operator_branch": "EMPIRICAL_SHARED_AFFINE",
        "source_sha256": {"DNA01_OBJECT_TABLE.tsv": sha256_file(args.object_table), "SOURCE_FASTA": sha256_file(args.fasta)},
        "runtime": {"python": sys.version.split()[0], "torch": torch.__version__, "numpy": np.__version__},
        "attention_backend": "standard_fallback" if args.model == "dnabert2" else "model_default",
        "repeat_policy": "bitwise exact required" if args.strict_repeat else "reported",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{args.output.name}.staging.", dir=args.output.parent))
    try:
        np.savez_compressed(stage / "EMBEDDINGS.npz", object_ids=np.asarray(object_ids), vectors=primary)
        np.savez_compressed(stage / "REPEAT_EMBEDDINGS.npz", object_ids=np.asarray(object_ids), vectors=repeat)
        (stage / "MODEL_INTERFACE.json").write_text(json.dumps(interface, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        validation = {"status": "PASS", "objects": len(object_ids), "dimension": int(primary.shape[1]), "windows": window_count, "repeat_windows": repeat_window_count, "repeat_bitwise_identical": repeat_exact, "repeat_max_abs_difference": max_delta, "all_finite": True, "elapsed_seconds": round(time.time() - started, 3)}
        (stage / "FORWARD_VALIDATION.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        write_metadata(stage, ["EMBEDDINGS.npz", "REPEAT_EMBEDDINGS.npz", "MODEL_INTERFACE.json", "FORWARD_VALIDATION.json"])
        os.replace(stage, args.output)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps(validation, sort_keys=True))


if __name__ == "__main__":
    main()
