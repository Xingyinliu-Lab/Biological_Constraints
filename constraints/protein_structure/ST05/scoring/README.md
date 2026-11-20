# ST05 scoring

## Environment

Create the standalone scoring environment described in [`ENVIRONMENT.md`](ENVIRONMENT.md). It is intentionally separate from model CUDA environments.

## Programs

### `score_st05.py`

Score any compatible invariant ST05 representation with the frozen numeric law.

```bash
python constraints/protein_structure/ST05/scoring/score_st05.py \
  --relations <RELATIONS> \
  --refset-release <REFSET_RELEASE> \
  --embeddings <EMBEDDINGS> \
  --model-id <MODEL_ID> \
  --output-dtype <OUTPUT_DTYPE> \
  --out <OUT>
```

### `st05_f_only_core.py`

Frozen numeric core imported by `score_st05.py`; it exposes no independent CLI.

### `st05_scoring_core.py`

Internal representation validation and invariant-law support imported by
`score_st05.py`; it is not a second scorer and exposes no CLI.

## Inputs

Use the exact RefSet paths and representation schema documented in [`../README.md`](../README.md). Inputs must cover the complete formal ID universe, contain finite numeric values and preserve frozen folds/dependency units.

## Outputs

Write outputs to `work/<constraint>/<model>/` for a new run. The final result
contains the primary score, uncertainty interval, population/coverage counts,
numerical settings, input identities and base/parent/cluster contributions.
`ST05_TRANSFORM_RESULT_LEDGER.tsv` is a required formal audit output: each
relation's residual, scale, tolerance, margin, `F_credit` and threshold rule
directly support the reported F hierarchy.

## One public scoring pipeline
