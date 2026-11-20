# PS04 × esmc-600m: forward implementation

This directory contains the actual model-by-constraint input construction, forward, coordinate mapping and representation-export implementation used for the registered result. It does not define or implement the formal score.

## Authorities

- Constraint, RefSet and representation contract: [`constraints/protein_sequence/PS04/README.md`](../../../../../../constraints/protein_sequence/PS04/README.md)
- Model checkpoint and environment: [`paper_cases/model_resources/protein_sequence/esmc-600m/README.md`](../../../../../model_resources/protein_sequence/esmc-600m/README.md)
- Unique public scoring handoff: [`scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)

## Forward scripts

- [`ps04_forward.py`](scripts/ps04_forward.py)
- [`run_forward.sh`](scripts/run_forward.sh)
- [`run_ps04_safe_cosine_preflight.py`](scripts/run_ps04_safe_cosine_preflight.py)

## Current portable invocation

```bash
# Run from the repository root.
python paper_cases/protein_sequence/PS04/models/esmc-600m/forward/scripts/ps04_forward.py esmc --fasta constraints/protein_sequence/PS04/refset/FORMAL_ENDPOINTS.faa --checkpoint /models/esmc-600m --workdir work/PS04/esmc-600m/forward
python constraints/protein_sequence/PS04/scoring/score_ps04.py --representations work/PS04/esmc-600m/forward/esmc/embeddings.safetensors --representation-index work/PS04/esmc-600m/forward/esmc/embedding_index.tsv --model-id esmc-600m --output-dir work/PS04/esmc-600m/score
```

## RunPlan and executed configuration

- The scripts above contain the released model parameters.

## Scoring boundary

Forward output must satisfy the public representation contract documented by the Constraint. The same public scorer is used for registered paper models and third-party models.
