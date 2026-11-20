"""Frozen JanusDNA checkpoint loader shared by paper Forward cases."""
from __future__ import annotations

import sys
import types
from pathlib import Path

import torch
import yaml

MODEL_LENGTH = 1024
BASE_TO_ID = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}
PAD_ID = 4


def load_model(source: Path, config_path: Path, checkpoint_path: Path):
    sys.path.insert(0, str(source))
    sys.modules.setdefault("wandb", types.SimpleNamespace())
    from janusdna.configuration_janusdna import JanusDNAConfig
    from janusdna.modeling_janusdna import JanusDNAForCausalLM

    JanusDNAForCausalLM._supports_flex_attn = True
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))["model"]["config"]
    config.pop("_target_", None)
    config.update(attn_implementation="flex_attention", output_hidden_states=False, use_mamba_kernels=False, attention_mask=False)
    model = JanusDNAForCausalLM(JanusDNAConfig(**config))
    state = torch.load(checkpoint_path, map_location="cpu", weights_only=False)["state_dict"]
    torch.nn.modules.utils.consume_prefix_in_state_dict_if_present(state, "model.")
    remapped = {}
    for direction in ("fwd", "bwd"):
        legacy = f"model.layers.4.attn.self_attn_{direction}.o_proj.weight"
        current = f"model.layers.4.attn.self_attn_{direction}.o_projs.0.weight"
        if legacy not in state or current in state or tuple(state[legacy].shape) != tuple(model.state_dict()[current].shape):
            raise RuntimeError(f"unexpected JanusDNA checkpoint projection keys: {legacy} -> {current}")
        state[current] = state.pop(legacy)
        remapped[legacy] = current
    missing, unexpected = model.load_state_dict(state, strict=False)
    ignored_metadata = sorted(key for key in unexpected if "torchmetrics" in key)
    unexpected_parameters = sorted(key for key in unexpected if key not in ignored_metadata)
    if missing or unexpected_parameters:
        raise RuntimeError(f"JanusDNA checkpoint mismatch after compatibility remap: missing={sorted(missing)} unexpected={unexpected_parameters}")
    model.eval().cuda().half()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model, sorted(missing), unexpected_parameters, remapped, ignored_metadata
