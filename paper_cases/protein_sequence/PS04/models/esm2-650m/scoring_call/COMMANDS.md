# Public scoring call

```bash
python constraints/protein_sequence/PS04/scoring/score_ps04.py \
  --representations work/PS04/esm2-650m/forward/esm2/embeddings.safetensors \
  --representation-index work/PS04/esm2-650m/forward/esm2/embedding_index.tsv \
  --model-id esm2-650m \
  --output-dir work/PS04/esm2-650m/score
```
