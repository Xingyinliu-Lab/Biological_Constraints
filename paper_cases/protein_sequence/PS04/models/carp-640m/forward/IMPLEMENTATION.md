# PS04 × carp-640m: forward implementation

This directory contains the actual model-by-constraint input construction, forward, coordinate mapping and representation-export implementation used for the registered result. It does not define or implement the formal score.

## Authorities

- Constraint, RefSet and representation contract: [`constraints/protein_sequence/PS04/README.md`](../../../../../../constraints/protein_sequence/PS04/README.md)
- Model checkpoint and environment: [`paper_cases/model_resources/protein_sequence/carp-640m/README.md`](../../../../../model_resources/protein_sequence/carp-640m/README.md)
- Unique public scoring handoff: [`scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)

## Forward scripts

- [`build_ps04_carp_anchor_trace.py`](scripts/build_ps04_carp_anchor_trace.py)
- [`ps04_carp640_forward.py`](scripts/ps04_carp640_forward.py)
- [`run.sh`](scripts/run.sh)
- [`run_ps04_safe_cosine_preflight.py`](scripts/run_ps04_safe_cosine_preflight.py)

## Current portable invocation

Provide the standard `TORCH_HOME/hub/checkpoints/carp_640M.pt` checkpoint path; the script derives `TORCH_HOME` from it for the CARP loader.

```bash
# Run from the repository root.
python paper_cases/protein_sequence/PS04/models/carp-640m/forward/scripts/ps04_carp640_forward.py carp --fasta constraints/protein_sequence/PS04/refset/FORMAL_ENDPOINTS.faa --checkpoint /torch-cache/hub/checkpoints/carp_640M.pt --workdir work/PS04/carp-640m/forward
python constraints/protein_sequence/PS04/scoring/score_ps04.py --representations work/PS04/carp-640m/forward/carp/embeddings.safetensors --representation-index work/PS04/carp-640m/forward/carp/embedding_index.tsv --model-id carp-640m --output-dir work/PS04/carp-640m/score
```

## RunPlan and executed configuration

- The scripts above contain the released model parameters.

## Scoring boundary

Forward output must satisfy the public representation contract documented by the Constraint. The same public scorer is used for registered paper models and third-party models.
