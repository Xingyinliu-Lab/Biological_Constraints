# Public scoring call

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings work/DNA02-02/dnabert-2/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA02-02/dnabert-2/forward/REPEAT_EMBEDDINGS.npz \
  --model-id dnabert-2 \
  --output work/DNA02-02/dnabert-2/score
```

This case does not contain another scorer. The command calls the unique public Constraint implementation.
