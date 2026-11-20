
# Public DNA01 scorer call

This case uses the unique scorer owned by `constraints/dna/DNA01/`:

```bash
python constraints/dna/DNA01/scoring/score_dna01.py \
  --embeddings work/DNA01/caduceus-ps/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA01/caduceus-ps/forward/REPEAT_EMBEDDINGS.npz \
  --model-key caduceus-ps \
  --output work/DNA01/caduceus-ps/score
```

The checked-in `output/DNA01_v0.5.4/` is the deterministic F-only materialization of the complete retained measurement ledger. No model Forward or affine refit was performed for this interface migration, and no result-local scorer exists in this case.
