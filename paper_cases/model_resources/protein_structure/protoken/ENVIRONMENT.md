# protoken executed environment

This file documents the environment observed for the published forward. It is model-owned provenance, not a project-wide environment and not a newly promised upstream lockfile.

## Runtime identity

- Evidence status: `READ_ONLY_VERIFIED_ON_EXECUTION_HOST`
- Verified host label: `kt`
- Verification date: `2026-10-16`
- Executed Python path: `/home/mik/u8t/ai_genome/m4_encoder_audit/env_protoken/bin/python3.10`
- Python version: `3.10.21`
- CUDA reported by PyTorch: `None`

## Packages used or relevant to forward execution

| Package | Executed version |
|---|---|
| `biopython` | `1.84` |
| `dm-haiku` | `0.0.12` |
| `jax` | `0.4.28` |
| `jaxlib` | `0.4.28` |
| `numpy` | `2.2.1` |
| `pandas` | `2.2.3` |
| `scikit-learn` | `1.5.1` |
| `scipy` | `1.15.1` |

Imports observed in released execution code: `jax`, `numpy`.

## Recreate the environment

Use the upstream model project's installation instructions, then match the recorded Python/package/CUDA versions where required by the selected execution adapter. Obtain the exact checkpoint and loading method from [`checkpoint_identity.json`](checkpoint_identity.json).

The historical Python path above is evidence only. Do not reproduce the author directory layout, and do not merge this model environment with another model's environment merely for convenience.
