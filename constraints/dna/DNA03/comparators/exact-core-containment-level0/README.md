# DNA03 exact-core-containment Level-0 reference

This deterministic reference is owned by [DNA03](../../README.md). It evaluates
the same frozen core-local relation and aggregation hierarchy as the public
model scorer. It is a white-box ceiling, not a foundation-model adapter.

## Run

```bash
python constraints/dna/DNA03/comparators/exact-core-containment-level0/score_dna03_exact_containment_level0.py \
  --output work/DNA03/exact-core-containment-level0
```

The program uses the released DNA03 Forward-view FASTA and relation tables, applies the frozen 20,000-draw dependency-cluster bootstrap, and writes `INDIVIDUAL_SCORES.tsv`, query and cluster-transition ledgers, `FORMAL_RESULT.json`, and `OUTPUT_MANIFEST.json`. Canonical outputs are retained in `results/` in this directory.
