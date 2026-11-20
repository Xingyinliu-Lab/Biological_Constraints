#!/usr/bin/env python3
"""Development-only PS02A endpoint encoding with the three PS01 protein conventions."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import torch
from safetensors.torch import save_file

BASE = Path.cwd()
ROOT = None
SOURCE = None
MODELS = {
    'esm2': None,
    'prott5': None,
    'esmc': None,
}

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    seqs = {}
    with gzip.open(path, 'rt') as f:
        key = None
        chunks = []
        for line in f:
            if line.startswith('>'):
                if key in wanted:
                    seqs[key] = ''.join(chunks)
                key = line[1:].split()[0]
                chunks = []
            elif key in wanted:
                chunks.append(line.strip())
        if key in wanted:
            seqs[key] = ''.join(chunks)
    return seqs

def prepare() -> None:
    manifest = ROOT / 'DEVELOPMENT_ENCODING_INPUT_MANIFEST.tsv'
    rows = list(csv.DictReader(manifest.open(newline=''), delimiter='\t'))
    assert len(rows) == 13720 and len({r['occurrence_id'] for r in rows}) == len(rows)
    seqs = read_fasta(SOURCE, {r['occurrence_id'] for r in rows})
    failures = []
    for r in rows:
        ident = r['occurrence_id']
        seq = seqs.get(ident)
        if seq is None:
            failures.append((ident, 'SOURCE_SEQUENCE_MISSING'))
        elif hashlib.sha256(seq.encode('ascii')).hexdigest() != r['sequence_sha256']:
            failures.append((ident, 'SOURCE_SEQUENCE_HASH_MISMATCH'))
    lengths = [len(s) for s in seqs.values()]
    report = {
        'status': 'READY' if not failures else 'BLOCKED',
        'kind': 'DEVELOPMENT_MODEL_FORWARD_INPUT',
        'manifest_sha256': sha(manifest), 'source_fasta_sha256': sha(SOURCE),
        'expected': len(rows), 'found': len(seqs), 'failures': failures[:15],
        'min_length': min(lengths) if lengths else None,
        'max_length': max(lengths) if lengths else None,
        'over_1022': sum(x > 1022 for x in lengths),
    }
    (ROOT / 'PREPARE_REPORT.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))
    if failures:
        raise SystemExit(1)
    with (ROOT / 'endpoint_sequences.faa').open('w') as out:
        for r in rows:
            ident = r['occurrence_id']
            out.write(f'>{ident}\n{seqs[ident]}\n')

def sequences() -> dict[str, str]:
    result = {}
    key = None
    with (ROOT / 'endpoint_sequences.faa').open() as f:
        for line in f:
            if line.startswith('>'):
                key = line[1:].strip()
                result[key] = ''
            else:
                result[key] += line.strip()
    return result

def batches(seqs: dict[str, str], max_count: int, max_tokens: int):
    batch, total = [], 0
    for ident in sorted(seqs, key=lambda k: (len(seqs[k]), k)):
        cost = len(seqs[ident]) + 2
        if batch and (len(batch) >= max_count or total + cost > max_tokens):
            yield batch
            batch, total = [], 0
        batch.append(ident)
        total += cost
    if batch:
        yield batch

def encode(model_name: str) -> None:
    seqs = sequences()
    outdir = ROOT / model_name
    outdir.mkdir(exist_ok=True)
    output = outdir / 'embeddings.safetensors'
    if output.exists():
        raise SystemExit(f'Existing output preserved: {output}')
    torch.cuda.set_device(0)
    torch.cuda.reset_peak_memory_stats()
    torch.manual_seed(20260929)
    checkpoint = MODELS[model_name]
    if model_name == 'esm2':
        sys.path.insert(0, str(BASE / 'm3_smoke/python_deps'))
        import esm
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(checkpoint))
        model = model.eval().half().cuda()
        convert = alphabet.get_batch_converter()
        def forward(ids):
            _, _, token = convert([(i, seqs[i]) for i in ids])
            with torch.inference_mode():
                h = model(token.cuda(), repr_layers=[33], return_contacts=False)['representations'][33]
            return [h[j, 1:len(seqs[i])+1].float().mean(0).cpu().contiguous() for j, i in enumerate(ids)]
        batch_iter = batches(seqs, 16, 4096)
        metadata = {'model': 'esm2_t33_650M_UR50D', 'layer': '33', 'pooling': 'residue_mean'}
    elif model_name == 'prott5':
        sys.path.insert(0, str(BASE / 'm3_encoder_audit/python_deps_prott5_444'))
        from transformers import T5EncoderModel, T5Tokenizer
        tok = T5Tokenizer.from_pretrained(str(checkpoint), do_lower_case=False, legacy=True, local_files_only=True)
        model = T5EncoderModel.from_pretrained(str(checkpoint), torch_dtype=torch.float16, local_files_only=True).cuda().eval()
        def forward(ids):
            spaced = [' '.join(re.sub(r'[UZOB]', 'X', seqs[i])) for i in ids]
            encoded = tok(spaced, add_special_tokens=True, padding=True, return_tensors='pt')
            with torch.inference_mode(), torch.autocast('cuda', dtype=torch.float16):
                h = model(input_ids=encoded['input_ids'].cuda(), attention_mask=encoded['attention_mask'].cuda()).last_hidden_state
            return [h[j, :len(seqs[i])].float().mean(0).cpu().contiguous() for j, i in enumerate(ids)]
        batch_iter = batches(seqs, 8, 2048)
        metadata = {'model': 'Rostlab/prot_t5_xl_half_uniref50-enc', 'pooling': 'residue_mean', 'UO_mapping': 'X'}
    else:
        from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer
        tok = EsmcTokenizer()
        model = EsmcForMaskedLM.from_pretrained(str(checkpoint), local_files_only=True).cuda().eval()
        def forward(ids):
            encoded = tok([seqs[i] for i in ids], return_tensors='pt', padding=True)
            encoded = {k: v.cuda() for k, v in encoded.items()}
            with torch.inference_mode(), torch.autocast('cuda', dtype=torch.float16):
                h = model(**encoded).last_hidden_state
            return [h[j, 1:len(seqs[i])+1].float().mean(0).cpu().contiguous() for j, i in enumerate(ids)]
        batch_iter = batches(seqs, 16, 4096)
        metadata = {'model': 'biohub/ESMC-600M', 'source_commit': 'bf343ba264b650dff7a073643725f9aaa1fdbe8d', 'pooling': 'residue_mean'}
    vectors, idx, failed = {}, [], []
    started = time.time()
    for number, ids in enumerate(batch_iter, 1):
        try:
            result = forward(ids)
        except Exception as exc:
            if len(ids) > 1:
                # Preserve every failed batch member and retry individually.
                for ident in ids:
                    try:
                        one = forward([ident])[0]
                        key = f'embedding_{len(vectors)+1:06d}'
                        vectors[key] = one
                        idx.append((ident, key, len(seqs[ident]), one.numel(), hashlib.sha256(one.numpy().tobytes()).hexdigest()))
                    except Exception as one_exc:
                        failed.append((ident, type(one_exc).__name__, str(one_exc)[:180]))
                continue
            failed.append((ids[0], type(exc).__name__, str(exc)[:180]))
            continue
        for ident, vec in zip(ids, result):
            if not torch.isfinite(vec).all():
                failed.append((ident, 'NONFINITE_EMBEDDING', ''))
                continue
            key = f'embedding_{len(vectors)+1:06d}'
            vectors[key] = vec
            idx.append((ident, key, len(seqs[ident]), vec.numel(), hashlib.sha256(vec.numpy().tobytes()).hexdigest()))
        if number % 20 == 0:
            print(json.dumps({'model': model_name, 'batches_done': number, 'encoded': len(idx), 'failed': len(failed)}), flush=True)
    if vectors:
        save_file(vectors, str(output), metadata=metadata)
    with (outdir / 'embedding_index.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(('occurrence_id', 'embedding_key', 'sequence_length', 'embedding_dim', 'embedding_sha256'))
        w.writerows(idx)
    with (outdir / 'failure_ledger.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(('occurrence_id', 'failure_type', 'detail'))
        w.writerows(failed)
    report = {'status': 'PASS' if len(idx) == len(seqs) and not failed else 'PARTIAL',
              'kind': 'DEVELOPMENT_MODEL_FORWARD', 'model': model_name,
              'endpoint_count': len(seqs), 'encoded_count': len(idx), 'failure_count': len(failed),
              'checkpoint_sha256': sha(checkpoint) if checkpoint.is_file() else sha(checkpoint/'pytorch_model.bin') if model_name == 'prott5' else sha(checkpoint/'model.safetensors'),
              'embedding_file_sha256': sha(output) if output.exists() else None,
              'embedding_dim': idx[0][3] if idx else None,
              'seconds': time.time()-started, 'peak_gpu_memory_mib': torch.cuda.max_memory_allocated()/2**20,
              'formal_result': False, 'formal_model_interface_release': False}
    (outdir / 'FORWARD_REPORT.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('command', choices=['prepare', 'esm2', 'prott5', 'esmc'])
    p.add_argument('--fasta', type=Path)
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--workdir', type=Path)
    args = p.parse_args()
    if any((args.fasta, args.checkpoint, args.workdir)):
        if args.command not in MODELS or not all((args.fasta, args.checkpoint, args.workdir)):
            p.error('current Forward requires a model command plus --fasta, --checkpoint and --workdir')
        ROOT = args.workdir
        ROOT.mkdir(parents=True, exist_ok=True)
        (ROOT / 'endpoint_sequences.faa').write_bytes(args.fasta.read_bytes())
        MODELS[args.command] = args.checkpoint
    if args.command == 'prepare': prepare()
    else: encode(args.command)
