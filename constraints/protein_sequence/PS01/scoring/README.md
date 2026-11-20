# PS01 scoring

## Environment

Create the standalone scoring environment described in [`ENVIRONMENT.md`](ENVIRONMENT.md). It is intentionally separate from model CUDA environments.

## Programs

### `score_ps01.py`

Evaluate compatible whole-sequence representations on the frozen PS01 RefSet.

```bash
python constraints/protein_sequence/PS01/scoring/score_ps01.py \
  --refset-dir <REFSET_DIR> \
  --representations <REPRESENTATIONS> \
  --model-id <MODEL_ID> \
  --output-dir <OUTPUT_DIR>
```

## Inputs

Use the exact RefSet paths and representation schema documented in [`../README.md`](../README.md). Inputs must cover the complete formal ID universe, contain finite numeric values and preserve frozen folds/dependency units.

## Outputs

Write outputs to `work/<constraint>/<model>/` for a new run. A final result contains F, its 95% interval, population/coverage counts, numerical settings, input identities, the required `INDIVIDUAL_SCORES.tsv.gz` F-credit ledger and `F_BOOTSTRAP.npy`. Temporary caches may be discarded; formal ledgers and bootstrap evidence must be retained with the result.

## One public scoring pipeline
