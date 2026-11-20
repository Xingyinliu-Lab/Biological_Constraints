# Public scoring call

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings work/DNA02-02/crossdna-8m/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA02-02/crossdna-8m/forward/REPEAT_EMBEDDINGS.npz \
  --model-id crossdna-8m \
  --output work/DNA02-02/crossdna-8m/score
```

This case does not contain another scorer. The command calls the unique public Constraint implementation.
