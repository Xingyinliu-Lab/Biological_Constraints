# RNA05 public scoring

`score_rna05.py` is the single authoritative model-scoring implementation. It accepts nucleotide-addressable model output and performs every frozen RNA05 operation after Forward.

```bash
python constraints/rna/RNA05/scoring/score_rna05.py \
  --vectors work/RNA05/my_model/REPRESENTATIONS.npz \
  --model-id my_model \
  --output work/RNA05/my_model/score
```

The optional `--ledger`, `--blocks` and `--fasta` inputs default to the complete released population and must be supplied together for a compatible alternative RefSet. The scorer validates their biological relationships, sequence/coordinate identities and numerical input types before fitting. The 20 outer leave-one-cluster-out fits, four inner folds, lambda grid, hierarchical weights, aggregation, seed and 20,000 bootstrap draws are fixed and cannot be overridden in a formal run.

Input schemas and every output are documented in the [Constraint README](../README.md). Environment requirements are in [`ENVIRONMENT.md`](ENVIRONMENT.md).

The scorer validates the frozen biological `(sequence_sha256, position0)` targets. Model-specific token-to-nucleotide mapping is performed and documented by the user during Forward; it is not implemented by the scorer.
