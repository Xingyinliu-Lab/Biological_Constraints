# carp-640m executed environment

This file documents the environment observed for the published forward. It is model-owned provenance, not a project-wide environment and not a newly promised upstream lockfile.

## Runtime identity

- Evidence status: `READ_ONLY_VERIFIED_ON_EXECUTION_HOST`
- Verified host label: `kt`
- Verification date: `2026-10-16`
- Executed Python path: `/home/mik/u8t/ai_genome/third_models_20261001/carp_640m_v0.1.0/env/bin/python`
- Python version: `3.10.20`
- CUDA reported by PyTorch: `13.0`

## Packages used or relevant to forward execution

| Package | Executed version |
|---|---|
| `biopython` | `1.85` |
| `biotite` | `0.41.2` |
| `huggingface-hub` | `0.36.2` |
| `jsonschema` | `4.26.0` |
| `numpy` | `1.26.4` |
| `pandas` | `2.3.3` |
| `safetensors` | `0.8.0` |
| `scikit-learn` | `1.7.2` |
| `scipy` | `1.15.3` |
| `torch` | `2.13.0` |
| `transformers` | `4.57.6` |

Imports observed in released execution code: `esm`, `numpy`, `safetensors`, `sequence_models`, `sklearn`, `torch`, `transformers`.

## Recreate the environment

Use the upstream model project's installation instructions, then match the recorded Python/package/CUDA versions where required by the selected execution adapter. Obtain the exact checkpoint and loading method from [`checkpoint_identity.json`](checkpoint_identity.json).

The historical Python path above is evidence only. Do not reproduce the author directory layout, and do not merge this model environment with another model's environment merely for convenience.
