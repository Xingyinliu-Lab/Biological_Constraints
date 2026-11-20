# prott5 executed environment

This file documents the environment observed for the published forward. It is model-owned provenance, not a project-wide environment and not a newly promised upstream lockfile.

## Runtime identity

- Evidence status: `READ_ONLY_VERIFIED_ON_EXECUTION_HOST`
- Verified host label: `kt`
- Verification date: `2026-10-16`
- Executed Python path: `/home/mik/u8t/ai_genome/third_models_20260914/env_esmc_py312/bin/python`
- Python version: `3.12.13`
- CUDA reported by PyTorch: `13.0`

## Packages used or relevant to forward execution

| Package | Executed version |
|---|---|
| `biopython` | `1.88` |
| `biotite` | `1.7.1` |
| `huggingface-hub` | `0.36.2` |
| `numpy` | `2.5.3` |
| `pandas` | `3.0.5` |
| `safetensors` | `0.8.0` |
| `scikit-learn` | `1.9.1` |
| `scipy` | `1.18.1` |
| `torch` | `2.11.0` |
| `transformers` | `4.57.6` |

Imports observed in released execution code: `esm`, `numpy`, `safetensors`, `sequence_models`, `sklearn`, `torch`, `transformers`.

## Recreate the environment

Use the upstream model project's installation instructions, then match the recorded Python/package/CUDA versions where required by the selected execution adapter. Obtain the exact checkpoint and loading method from [`checkpoint_identity.json`](checkpoint_identity.json).

The historical Python path above is evidence only. Do not reproduce the author directory layout, and do not merge this model environment with another model's environment merely for convenience.
