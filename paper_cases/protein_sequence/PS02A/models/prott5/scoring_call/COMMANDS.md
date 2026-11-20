# Public scoring call

```bash
python constraints/protein_sequence/PS02A/scoring/score_ps02a.py \
  --representations work/PS02A/prott5/forward/prott5/embeddings.safetensors \
  --representation-index work/PS02A/prott5/forward/prott5/embedding_index.tsv \
  --model-id prott5 \
  --output-dir work/PS02A/prott5/score
```

The scorer fixes all PS02A scientific parameters and writes raw measurements,
ID-level F credits, F and F 95% CI, coverage and output identities.
