# ST02 scoring

## Environment

Create the standalone scoring environment described in [`ENVIRONMENT.md`](ENVIRONMENT.md). It is intentionally separate from model CUDA environments.

## Programs

### `score_st02.py`

Evaluate compatible structure-model embeddings with the frozen ST02 method.

```bash
python constraints/protein_structure/ST02/scoring/score_st02.py \
  --refset-dir <REFSET_DIR> \
  --embeddings <EMBEDDINGS> \
  --model-id <MODEL_ID> \
  --output-dir <OUTPUT_DIR>
```

### `st02_numeric_core.py`

Frozen contrast, hierarchy and bootstrap implementation imported by
`score_st02.py`; it is not a second scorer and exposes no CLI.

## Inputs

Use the exact RefSet paths and representation schema documented in [`../README.md`](../README.md). Inputs must cover the complete formal ID universe, contain finite numeric values and preserve frozen folds/dependency units.

## Outputs

Write outputs to `work/<constraint>/<model>/score/` for a new run. The public
output consists of `SCORE_SUMMARY.json`, `INDIVIDUAL_SCORES.tsv`, hierarchy F
ledgers, `F_BOOTSTRAP.npy`, and `OUTPUT_MANIFEST.json`. ID-level measurements
and F scores are formal public outputs, not disposable intermediates.

## One public scoring pipeline
