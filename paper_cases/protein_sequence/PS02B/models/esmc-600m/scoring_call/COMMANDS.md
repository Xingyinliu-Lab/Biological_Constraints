# Public scoring call

```bash
python constraints/protein_sequence/PS02B/scoring/score_ps02b.py \
  --representations work/PS02B/esmc-600m/forward/embeddings.safetensors \
  --representation-index work/PS02B/esmc-600m/forward/embedding_index.tsv \
  --model-id esmc-600m \
  --output-dir work/PS02B/esmc-600m/score
```

The Constraint fixes all scoring parameters and writes raw measurements, ID-level F credits, formal F/F 95% CI, coverage and output identities.
