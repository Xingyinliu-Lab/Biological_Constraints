# Public scoring call

```bash
python constraints/protein_sequence/PS02A/scoring/score_ps02a.py \
  --representations work/PS02A/esm3-seq/forward/embeddings.safetensors \
  --representation-index work/PS02A/esm3-seq/forward/embedding_index.tsv \
  --model-id esm3-seq \
  --output-dir work/PS02A/esm3-seq/score
```

The scorer fixes all PS02A scientific parameters and writes raw measurements,
ID-level F credits, F and F 95% CI, coverage and output identities.
