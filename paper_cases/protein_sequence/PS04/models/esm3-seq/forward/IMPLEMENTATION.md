# PS04 × esm3-seq: forward implementation

This directory contains the actual model-by-constraint input construction, forward, coordinate mapping and representation-export implementation used for the registered result. It does not define or implement the formal score.

## Authorities

- Constraint, RefSet and representation contract: [`constraints/protein_sequence/PS04/README.md`](../../../../../../constraints/protein_sequence/PS04/README.md)
- Model checkpoint and environment: [`paper_cases/model_resources/protein_sequence/esm3-seq/README.md`](../../../../../model_resources/protein_sequence/esm3-seq/README.md)
- Unique public scoring handoff: [`scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)

## Forward scripts

- [`esm3_seq_mean_extract.py`](scripts/esm3_seq_mean_extract.py)

## Current portable invocation

Run from the repository root. `endpoint_id` is the first column written to the representation index and preserves the frozen FASTA object identity.

```bash
python paper_cases/protein_sequence/PS04/models/esm3-seq/forward/scripts/esm3_seq_mean_extract.py --fasta constraints/protein_sequence/PS04/refset/FORMAL_ENDPOINTS.faa --id-column endpoint_id --output-dir work/PS04/esm3-seq/forward/esm3
python constraints/protein_sequence/PS04/scoring/score_ps04.py --representations work/PS04/esm3-seq/forward/esm3/embeddings.safetensors --representation-index work/PS04/esm3-seq/forward/esm3/embedding_index.tsv --model-id esm3-seq --output-dir work/PS04/esm3-seq/score
```

## RunPlan and executed configuration

- The scripts above contain the released model parameters.

## Scoring boundary

Forward output must satisfy the public representation contract documented by the Constraint. The same public scorer is used for registered paper models and third-party models.
