# Public scoring call

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings work/DNA02-02/caduceus-ps/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA02-02/caduceus-ps/forward/REPEAT_EMBEDDINGS.npz \
  --model-id caduceus-ps \
  --output work/DNA02-02/caduceus-ps/score
```

This case does not contain another scorer. The command calls the unique public Constraint implementation.
