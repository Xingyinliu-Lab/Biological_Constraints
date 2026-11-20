# Public scoring call

```bash
python constraints/protein_sequence/PS02B/scoring/score_ps02b.py \
  --representations work/PS02B/esm3-seq/forward/embeddings.safetensors \
  --representation-index work/PS02B/esm3-seq/forward/embedding_index.tsv \
  --model-id esm3-seq \
  --output-dir work/PS02B/esm3-seq/score
```

The Constraint fixes all scoring parameters and writes raw measurements, ID-level F credits, formal F/F 95% CI, coverage and output identities.
