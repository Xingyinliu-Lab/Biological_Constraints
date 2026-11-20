# rinalmo-giga executed environment

This file documents the environment observed for the published forward. It is model-owned provenance, not a project-wide environment and not a newly promised upstream lockfile.

## Runtime identity

- Evidence status: `READ_ONLY_VERIFIED_ON_EXECUTION_HOST`
- Verified host label: `kt`
- Verification date: `2026-10-16`
- Executed Python path: `/home/mik/fd/m2_rna/jobs/rna09_closure_20261005/rinalmo_official/venv/bin/python`
- Python version: `3.11.16`
- CUDA reported by PyTorch: `12.1`

## Packages used or relevant to forward execution

| Package | Executed version |
|---|---|
| `numpy` | `1.26.4` |
| `torch` | `2.1.2` |

Imports observed in released execution code: `numpy`, `rinalmo`, `torch`.

## Recreate the environment

Use the upstream model project's installation instructions, then match the recorded Python/package/CUDA versions where required by the selected execution adapter. Obtain the exact checkpoint and loading method from [`checkpoint_identity.json`](checkpoint_identity.json).

The historical Python path above is evidence only. Do not reproduce the author directory layout, and do not merge this model environment with another model's environment merely for convenience.
