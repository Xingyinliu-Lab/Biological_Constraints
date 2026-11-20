#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import os
import random
import shutil
import sys
import tempfile
import time
import types
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import yaml

WINDOW_BP = 1000
MODEL_LENGTH = 1024
BATCH_SIZE = 128
BASE_TO_ID = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}
PAD_ID = 4


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path):
    records = {}
    key = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            if key in records:
                raise ValueError("duplicate FASTA id")
            records[key] = ""
        else:
            if key is None:
                raise ValueError("sequence before FASTA header")
            records[key] += line.strip().upper()
    return records


def load_model(source, config_path, checkpoint_path):
    sys.path.insert(0, str(source))
    sys.modules.setdefault("wandb", types.SimpleNamespace())
    from janusdna.configuration_janusdna import JanusDNAConfig
    from janusdna.modeling_janusdna import JanusDNAForCausalLM

    JanusDNAForCausalLM._supports_flex_attn = True
    cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))["model"]["config"]
    cfg.pop("_target_", None)
    cfg.update(
        attn_implementation="flex_attention",
        output_hidden_states=False,
        use_mamba_kernels=False,
        attention_mask=False,
    )
    model = JanusDNAForCausalLM(JanusDNAConfig(**cfg))
    state = torch.load(checkpoint_path, map_location="cpu", weights_only=False)["state_dict"]
    torch.nn.modules.utils.consume_prefix_in_state_dict_if_present(state, "model.")
    remapped = {}
    for direction in ("fwd", "bwd"):
        legacy = f"model.layers.4.attn.self_attn_{direction}.o_proj.weight"
        current = f"model.layers.4.attn.self_attn_{direction}.o_projs.0.weight"
        if legacy not in state or current in state:
            raise RuntimeError(f"unexpected JanusDNA checkpoint projection keys: {legacy} -> {current}")
        if tuple(state[legacy].shape) != tuple(model.state_dict()[current].shape):
            raise RuntimeError(f"JanusDNA projection shape mismatch: {legacy} -> {current}")
        state[current] = state.pop(legacy)
        remapped[legacy] = current
    missing, unexpected = model.load_state_dict(state, strict=False)
    ignored_metadata = sorted(k for k in unexpected if "torchmetrics" in k)
    unexpected_parameters = sorted(k for k in unexpected if k not in ignored_metadata)
    if missing or unexpected_parameters:
        raise RuntimeError(
            f"JanusDNA checkpoint mismatch after compatibility remap: "
            f"missing={sorted(missing)} unexpected={unexpected_parameters}"
        )
    model.eval().cuda().half()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model, sorted(missing), unexpected_parameters, remapped, ignored_metadata


def encode(model, object_ids, sequences, run_name):
    tasks = []
    weights = defaultdict(int)
    sums = {}
    for object_id in object_ids:
        sequence = sequences[object_id]
        for start in range(0, len(sequence), WINDOW_BP):
            window = sequence[start : start + WINDOW_BP]
            tasks.append((object_id, window))
            weights[object_id] += len(window)

    for start in range(0, len(tasks), BATCH_SIZE):
        batch = tasks[start : start + BATCH_SIZE]
        encoded = []
        lengths = []
        for _, sequence in batch:
            try:
                ids = [BASE_TO_ID[base] for base in sequence]
            except KeyError as exc:
                raise ValueError(f"unsupported DNA character: {exc.args[0]}") from exc
            lengths.append(len(ids))
            encoded.append(ids + [PAD_ID] * (MODEL_LENGTH - len(ids)))
        while len(encoded) < BATCH_SIZE:
            encoded.append([PAD_ID] * MODEL_LENGTH)

        input_ids = torch.tensor(encoded, dtype=torch.long, device="cuda")
        with torch.inference_mode(), torch.autocast(device_type="cuda", dtype=torch.float16):
            hidden = model.model(
                input_ids=input_ids,
                output_hidden_states=False,
                return_dict=True,
                use_cache=False,
            ).last_hidden_state.float()
        for row, ((object_id, _), length) in enumerate(zip(batch, lengths)):
            vector = hidden[row, :length].mean(dim=0).cpu()
            sums[object_id] = sums.get(object_id, 0) + vector * length
        if start == 0 or start + len(batch) == len(tasks) or (start // BATCH_SIZE) % 50 == 0:
            print(
                json.dumps(
                    {
                        "model": "janusdna",
                        "pass": run_name,
                        "windows_done": min(start + len(batch), len(tasks)),
                        "windows_total": len(tasks),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
    result = torch.stack([(sums[oid] / weights[oid]).float() for oid in object_ids]).numpy()
    if not np.isfinite(result).all():
        raise ValueError("non-finite embeddings")
    return result, len(tasks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--object-table", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--single-pass", choices=("primary", "repeat"))
    parser.add_argument("--reload-between-passes", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("output path must not already exist")

    random.seed(20261014)
    np.random.seed(20261014)
    torch.manual_seed(20261014)
    torch.cuda.manual_seed_all(20261014)

    sequences = read_fasta(args.fasta)
    with args.object_table.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    object_ids = []
    for row in rows:
        object_id = row["object_id"]
        sequence = sequences.get(object_id)
        if (
            sequence is None
            or len(sequence) != int(row["length_bp"])
            or hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]
        ):
            raise ValueError("object identity mismatch")
        object_ids.append(object_id)
    if args.limit is not None:
        object_ids = object_ids[: args.limit]

    model, missing, unexpected, remapped, ignored_metadata = load_model(args.source, args.config, args.checkpoint)
    started = time.time()
    if args.single_pass:
        vectors, windows = encode(model, object_ids, sequences, args.single_pass)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        if args.output.exists():
            raise SystemExit("single-pass output must not already exist")
        stage_file = args.output.with_name(f".{args.output.name}.staging.{os.getpid()}")
        with stage_file.open("wb") as handle:
            np.savez_compressed(handle, object_ids=np.asarray(object_ids), vectors=vectors)
        os.replace(stage_file, args.output)
        print(json.dumps({"status": "PASS", "pass": args.single_pass, "objects": len(object_ids), "dimension": int(vectors.shape[1]), "windows": windows, "elapsed_seconds": round(time.time() - started, 3), "missing_keys": missing, "unexpected_parameter_keys": unexpected, "compatibility_key_remap": remapped, "ignored_checkpoint_metadata_keys": ignored_metadata}, sort_keys=True))
        return
    primary, primary_windows = encode(model, object_ids, sequences, "primary")
    if args.reload_between_passes:
        del model
        torch.cuda.empty_cache()
        model, missing, unexpected, remapped, ignored_metadata = load_model(args.source, args.config, args.checkpoint)
    repeat, repeat_windows = encode(model, object_ids, sequences, "repeat")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{args.output.name}.staging.", dir=args.output.parent))
    try:
        np.savez_compressed(stage / "EMBEDDINGS.npz", object_ids=np.asarray(object_ids), vectors=primary)
        np.savez_compressed(stage / "REPEAT_EMBEDDINGS.npz", object_ids=np.asarray(object_ids), vectors=repeat)
        interface = {
            "artifact": "DNA02_02_MODEL_INTERFACE",
            "status": "EXECUTED_FROZEN_INTERFACE",
            "model_key": "janusdna",
            "model_id": "JanusDNA_len1k_dmodel32_interdim128_8layers",
            "representation": "native JanusDNA final bidirectional-fusion hidden state; mean over real bases within each nonoverlapping 1000-bp window; bp-weighted window mean",
            "window_bp": WINDOW_BP,
            "model_input_length": MODEL_LENGTH,
            "tokenizer": "official CharacterTokenizer ids: A/C/G/T/N=7/8/9/10/11; PAD=4; no added special tokens",
            "precision": "model float16; pooling and released vectors float32",
            "operator_branch": "FROZEN_DNA02_02_FORWARD_INTERFACE",
            "checkpoint_sha256": sha(args.checkpoint),
            "source_sha256": {
                "DNA02_02_FORWARD_OBJECT_TABLE.tsv": sha(args.object_table),
                "SOURCE_FASTA": sha(args.fasta),
            },
            "runtime": {
                "python": sys.version.split()[0],
                "torch": torch.__version__,
                "numpy": np.__version__,
                "attention": "flex_attention",
                "use_mamba_kernels": False,
            },
            "checkpoint_load": {
                "status": "PASS_ALL_MODEL_PARAMETERS_LOADED",
                "missing_keys": missing,
                "unexpected_parameter_keys": unexpected,
                "compatibility_key_remap": remapped,
                "ignored_checkpoint_metadata_keys": ignored_metadata,
            },
        }
        (stage / "MODEL_INTERFACE.json").write_text(json.dumps(interface, indent=2, sort_keys=True) + "\n")
        repeat_delta = float(np.max(np.abs(primary - repeat)))
        validation = {
            "status": "PASS" if repeat_delta == 0.0 else "FAIL_REPEAT_MISMATCH",
            "objects": len(object_ids),
            "dimension": int(primary.shape[1]),
            "windows": primary_windows,
            "repeat_windows": repeat_windows,
            "checkpoint_reload_between_passes": bool(args.reload_between_passes),
            "repeat_max_abs_difference": repeat_delta,
            "all_finite": True,
            "elapsed_seconds": round(time.time() - started, 3),
        }
        if repeat_delta != 0.0:
            raise RuntimeError(json.dumps(validation, sort_keys=True))
        (stage / "FORWARD_VALIDATION.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
        names = ["EMBEDDINGS.npz", "REPEAT_EMBEDDINGS.npz", "MODEL_INTERFACE.json", "FORWARD_VALIDATION.json"]
        files = {name: {"sha256": sha(stage / name), "bytes": (stage / name).stat().st_size} for name in names}
        (stage / "OUTPUT_MANIFEST.json").write_text(json.dumps({"status": "PASS", "files": files}, indent=2, sort_keys=True) + "\n")
        names.append("OUTPUT_MANIFEST.json")
        (stage / "SHA256SUMS").write_text("".join(f"{sha(stage / name)}  {name}\n" for name in names))
        (stage / "RUN_COMPLETE").write_text("PASS\n")
        os.replace(stage, args.output)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps(validation, sort_keys=True))


if __name__ == "__main__":
    main()
