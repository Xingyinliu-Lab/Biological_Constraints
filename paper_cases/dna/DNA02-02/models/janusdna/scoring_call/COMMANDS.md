# Public scoring call

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings work/DNA02-02/janusdna/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA02-02/janusdna/forward/REPEAT_EMBEDDINGS.npz \
  --model-id janusdna \
  --output work/DNA02-02/janusdna/score
```

This case does not contain another scorer. The command calls the unique public Constraint implementation.
