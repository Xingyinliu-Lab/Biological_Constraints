# RNA04 public scoring

`score_rna04.py` is the single public F-only model-scoring implementation. It accepts model Forward output and completes all frozen RNA04 measurement, cross-fit, F aggregation and uncertainty calculations.

```bash
python constraints/rna/RNA04/scoring/score_rna04.py \
  --vectors work/RNA04/my_model/EMBEDDINGS.npz \
  --model-id my_model \
  --output work/RNA04/my_model/score
```

The optional `--ledger` and `--fasta` arguments default to the released population and must be supplied together for a compatible alternative RefSet. The scorer validates the biological relationships and records the actual input identities; it does not use the released files' byte hashes as a general admission gate. Users do not choose the fold count, regularization grid, hierarchy, bootstrap seed or number of draws; these are fixed by the release.

Inputs and every output column are described in the [Constraint README](../README.md). The scoring environment is documented in [`ENVIRONMENT.md`](ENVIRONMENT.md).
