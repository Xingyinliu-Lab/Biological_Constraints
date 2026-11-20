# dnabert-2 executed environment

This file documents the environment observed for the published forward. It is model-owned provenance, not a project-wide environment and not a newly promised upstream lockfile.

## Runtime identity

- Evidence status: `READ_ONLY_VERIFIED_ON_EXECUTION_HOST`
- Verified host label: `kt`
- Verification date: `2026-10-16`
- Executed Python path: `/home/mik/u8t/ai_genome/dna_models/envs_20261014/dnabert2_direct_hf/bin/python`
- Python version: `3.10.18`
- CUDA reported by PyTorch: `12.4`

## Packages used or relevant to forward execution

| Package | Executed version |
|---|---|
| `biotite` | `0.39.0` |
| `huggingface-hub` | `0.36.2` |
| `numpy` | `1.26.4` |
| `pandas` | `2.3.1` |
| `scipy` | `1.15.3` |
| `torch` | `2.5.1` |
| `transformers` | `4.29.2` |

Imports observed in released execution code: `janusdna`, `multimolecule`, `numpy`, `torch`, `transformers`, `typing_extensions`, `yaml`.

## Recreate the environment

Use the upstream model project's installation instructions, then match the recorded Python/package/CUDA versions where required by the selected execution adapter. Obtain the exact checkpoint and loading method from [`checkpoint_identity.json`](checkpoint_identity.json).

The historical Python path above is evidence only. Do not reproduce the author directory layout, and do not merge this model environment with another model's environment merely for convenience.
