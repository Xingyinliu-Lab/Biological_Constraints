# esm3-struct executed environment

This file documents the environment observed for the published forward. It is model-owned provenance, not a project-wide environment and not a newly promised upstream lockfile.

## Runtime identity

- Evidence status: `READ_ONLY_VERIFIED_ON_EXECUTION_HOST`
- Verified host label: `kt`
- Verification date: `2026-10-16`
- Executed Python path: `/media/datau28t/mik/ai_genome/envs/esm3-sm-open-v1/bin/python3.10`
- Python version: `3.10.22`
- CUDA reported by PyTorch: `12.4`

## Packages used or relevant to forward execution

| Package | Executed version |
|---|---|
| `biopython` | `1.88` |
| `biotite` | `0.41.2` |
| `huggingface-hub` | `0.36.2` |
| `numpy` | `1.26.4` |
| `pandas` | `2.3.3` |
| `safetensors` | `0.8.0` |
| `scikit-learn` | `1.7.2` |
| `scipy` | `1.15.3` |
| `torch` | `2.5.1+cu124` |
| `transformers` | `4.47.1` |

Imports observed in released execution code: `esm`, `numpy`, `torch`.

## Recreate the environment

Use the upstream model project's installation instructions, then match the recorded Python/package/CUDA versions where required by the selected execution adapter. Obtain the exact checkpoint and loading method from [`checkpoint_identity.json`](checkpoint_identity.json).

The historical Python path above is evidence only. Do not reproduce the author directory layout, and do not merge this model environment with another model's environment merely for convenience.
