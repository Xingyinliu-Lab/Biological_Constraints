# ST05 comparator: rigid-geometry white box

This deterministic reference evaluates the frozen ST05 rigid-transformation
relations from geometry. It is a white-box ceiling, not evidence of learned
structural knowledge.

## Run

```bash
python constraints/protein_structure/ST05/comparators/geometry-whitebox/score_st05_rigid_geometry_whitebox.py \
  --release constraints/protein_structure/ST05/refset \
  --out-dir work/ST05/geometry_whitebox
```

The implementation uses the fixed public float64 geometry tolerance and the
same F-only aggregation core as the model scorer.

## Registered outputs

- `results/WB_TRANSFORM_LEDGER.tsv`: ID-level geometric residuals, tolerances,
  normalized residuals and `F_credit` values for all 960 relations.
- `results/WB_RESULT.json`: registered white-box reference summary.
- `results/WB_VALIDATION.json`: coordinate, rotation, panel and finite-value QC.

The registered reference has `F=1`. Its inferential CI is not applicable
because it is a deterministic analytical ceiling.
