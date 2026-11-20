#!/usr/bin/env python3
"""DNA06 v0.1.1 extractor for NT-v2, DNABERT-2 and Caduceus-PS.

The model/environment loaders are inherited from the validated DNA01/DNA02
entry points.  Only the DNA03 readout is revised here.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, os, shutil, sys, tempfile, time
from pathlib import Path

import numpy as np
import torch

from dna06_continuous_mapping import (
    WindowProjection, frozen_windows, merge_window_projections,
    native_single_base_spans, pool_codon_local, project_tokens_to_nucleotides,
)
from dna06_forward_publication import publish_forward

MODEL_IDS = {
    "nt": "InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
    "dnabert2": "zhihan1996/DNABERT-2-117M",
    "caduceus_ps": "kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
}
MODEL_REVISIONS = {
    "nt": "f34324c6fde36a4f635f0f1f06cac5d25acd6798",
    "dnabert2": "7bce263b15377fc15361f52cfab88f8b586abda0",
    "caduceus_ps": "d89eeb853136ea64da7feb3d0c8e909771b17ae6",
}
BATCH = {"nt": 32, "dnabert2": 24, "caduceus_ps": 20}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path: Path):
    result, key = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]; result[key] = ""
        else:
            result[key] += line.strip().upper()
    return result


def read_rows(path: Path, fasta):
    with path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        seq = fasta[row["view_id"]]
        cs, ce = int(row["core_offset_start_0based"]), int(row["core_offset_end_exclusive"])
        if len(seq) != int(row["target_length_nt"]) or hashlib.sha256(seq.encode()).hexdigest() != row["view_sequence_sha256"] or ce - cs != 3:
            raise ValueError("view identity or core coordinates mismatch")
    return rows


def load_model(name):
    from transformers import AutoModel, AutoModelForMaskedLM, AutoTokenizer
    model_id = os.environ.get("DNA06_MODEL_PATH", MODEL_IDS[name])
    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True, local_files_only=True)
    if name == "dnabert2":
        model = AutoModel.from_pretrained(model_id, trust_remote_code=True, local_files_only=True)
        module = sys.modules[type(model).__module__]
        if hasattr(module, "flash_attn_qkvpacked_func"):
            module.flash_attn_qkvpacked_func = None
        model.float()
    else:
        model = AutoModelForMaskedLM.from_pretrained(model_id, trust_remote_code=True, local_files_only=True).half()
    model.eval().cuda()
    for parameter in model.parameters(): parameter.requires_grad_(False)
    return model_id, tokenizer, model


def build_tasks(rows, fasta):
    tasks = []
    for row in rows:
        view, seq = row["view_id"], fasta[row["view_id"]]
        for start, end in frozen_windows(len(seq)):
            tasks.append((end - start, view, start, end, seq[start:end]))
    return sorted(tasks, key=lambda x: (x[0], x[1], x[2]))


def tokenize(tokenizer, sequences, model_name):
    try:
        enc = tokenizer(sequences, return_tensors="pt", padding=True, truncation=False, return_offsets_mapping=True)
        offsets = enc.pop("offset_mapping").cpu().numpy()
        return enc, offsets
    except (NotImplementedError, TypeError, ValueError):
        if model_name not in {"nt", "caduceus_ps"}:
            raise ValueError("formal multi-base tokenizer must return authoritative offsets")
        enc = tokenizer(sequences, return_tensors="pt", padding=True, truncation=False)
        ids, special = enc["input_ids"], set(tokenizer.all_special_ids)
        offsets = np.zeros((ids.shape[0], ids.shape[1], 2), dtype=np.int64)
        for r, seq in enumerate(sequences):
            positions = [i for i, token in enumerate(ids[r].tolist()) if token not in special and int(enc.get("attention_mask", torch.ones_like(ids))[r, i])]
            if model_name == "caduceus_ps":
                if len(positions) != len(seq):
                    raise ValueError("native single-base token count mismatch")
                spans = native_single_base_spans(len(seq))
            else:
                cursor, spans, rebuilt = 0, [], []
                for i in positions:
                    token = tokenizer.convert_ids_to_tokens(int(ids[r, i]))
                    if not token or any(base not in "ACGT" for base in token):
                        raise ValueError("NT token cannot be mapped exactly to nucleotides")
                    spans.append((cursor, cursor + len(token))); cursor += len(token); rebuilt.append(token)
                if cursor != len(seq) or "".join(rebuilt) != seq:
                    raise ValueError("NT token-to-nucleotide reconstruction mismatch")
            for i, span in zip(positions, spans): offsets[r, i] = span
        return enc, offsets


def encode_pass(name, tokenizer, model, rows, fasta):
    tasks, by_view = build_tasks(rows, fasta), {}
    special = set(tokenizer.all_special_ids)
    for base in range(0, len(tasks), BATCH[name]):
        batch = tasks[base:base + BATCH[name]]
        enc, offsets = tokenize(tokenizer, [x[4] for x in batch], name)
        gpu = {k: v.cuda() for k, v in enc.items()}
        with torch.inference_mode(), torch.autocast(device_type="cuda", dtype=torch.float16, enabled=name != "dnabert2"):
            out = model(**gpu, output_hidden_states=name != "dnabert2", return_dict=True)
            hidden = (out.last_hidden_state if hasattr(out, "last_hidden_state") else out[0]) if name == "dnabert2" else out.hidden_states[-1]
        valid = gpu.get("attention_mask", torch.ones(hidden.shape[:2], device=hidden.device, dtype=torch.long)).bool()
        for token_id in special: valid &= gpu["input_ids"].ne(token_id)
        for i, (_, view, start, end, seq) in enumerate(batch):
            spans = []
            for ok, (a, b) in zip(valid[i].cpu().numpy().astype(bool), offsets[i]):
                spans.append((int(a), int(b)) if ok else (0, 0))
            positions, coverage = project_tokens_to_nucleotides(hidden[i].float().cpu().numpy(), spans, len(seq))
            by_view.setdefault(view, []).append(WindowProjection(start, end, positions, coverage))
    embeddings, keys, global_coverage = [], [], {}
    row_by_view = {r["view_id"]: r for r in rows}
    for view in [r["view_id"] for r in rows]:
        row, seq = row_by_view[view], fasta[view]
        merged, counts = merge_window_projections(len(seq), by_view[view])
        embedding = pool_codon_local(merged, int(row["core_offset_start_0based"]), int(row["core_offset_end_exclusive"]))
        keys.append(view); embeddings.append(embedding)
        global_coverage[view] = {"min": int(counts.min()), "max": int(counts.max())}
    return keys, np.stack(embeddings), len(tasks), global_coverage


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=sorted(MODEL_IDS), required=True)
    ap.add_argument("--membership", type=Path, required=True); ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True); ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    if args.output.exists(): raise SystemExit("output path exists")
    torch.use_deterministic_algorithms(True); torch.backends.cudnn.benchmark = False
    torch.backends.cuda.matmul.allow_tf32 = False; torch.backends.cudnn.allow_tf32 = False
    fasta = read_fasta(args.fasta); rows = read_rows(args.membership, fasta)
    if args.limit: rows = rows[:args.limit]
    model_id, tokenizer, model = load_model(args.model); started = time.time()
    ids, embeddings, task_count, coverage = encode_pass(args.model, tokenizer, model, rows, fasta)
    ids2, repeat_embeddings, repeat_tasks, coverage2 = encode_pass(args.model, tokenizer, model, rows, fasta)
    exact = ids == ids2 and np.array_equal(embeddings, repeat_embeddings) and coverage == coverage2
    if not exact: raise RuntimeError("strict repeat mismatch")
    publish_forward(args.output,model_key=args.model,model_identity={'model_id':MODEL_IDS[args.model],'model_revision':MODEL_REVISIONS[args.model],'model_path':model_id},
        source_paths={'membership':args.membership,'fasta':args.fasta},primary=(ids,embeddings),repeat=(ids2,repeat_embeddings),
        validation={'tasks':task_count,'repeat_tasks':repeat_tasks,'global_nucleotide_coverage':coverage,'repeat_global_nucleotide_coverage':coverage2,'elapsed_seconds':round(time.time()-started,3)},
        code_paths={'adapter':Path(__file__),'mapping':Path(__file__).with_name('dna06_continuous_mapping.py'),'publisher':Path(__file__).with_name('dna06_forward_publication.py')})


if __name__ == "__main__": main()
