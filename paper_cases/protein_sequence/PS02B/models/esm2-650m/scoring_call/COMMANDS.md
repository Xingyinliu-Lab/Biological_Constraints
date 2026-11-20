# Public scoring call

```bash
python constraints/protein_sequence/PS02B/scoring/score_ps02b.py \
  --representations work/PS02B/esm2-650m/forward/embeddings.safetensors \
  --representation-index work/PS02B/esm2-650m/forward/embedding_index.tsv \
  --model-id esm2-650m \
  --output-dir work/PS02B/esm2-650m/score
```

The Constraint fixes all scoring parameters and writes raw measurements, ID-level F credits, formal F/F 95% CI, coverage and output identities.
