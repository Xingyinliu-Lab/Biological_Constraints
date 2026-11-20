# DNA02-02 public scoring package

`score_dna0202.py` is the only public scoring entry. It accepts model Forward outputs and performs the complete frozen measurement chain: RefSet coverage validation, stratified repeat calibration, safe-cosine contrasts, ladder credits, hierarchical aggregation and outer-cluster bootstrap. `INDIVIDUAL_SCORES.tsv` retains all four component cosine measurements as well as their contrasts and credits.

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings <PRIMARY_NPZ> \
  --repeat-embeddings <REPEAT_NPZ> \
  --model-id <MODEL_ID> \
  --output <EMPTY_OUTPUT_DIRECTORY>
```

The scorer locates the released RefSet automatically. See [../README.md](../README.md) for the exact NPZ contract and output columns, and [ENVIRONMENT.md](ENVIRONMENT.md) for its independent CPU environment.

Internal modules have non-overlapping roles:

- `dna0202_numeric_core.py`: frozen mathematical primitives and statistics;
- `dna0202_io.py`: public input/output parsing and validation;
- `score_dna0202.py`: the single user-facing orchestration entry.

There is no model-specific or paper-only scorer. Paper cases and third-party models use this same entry.
