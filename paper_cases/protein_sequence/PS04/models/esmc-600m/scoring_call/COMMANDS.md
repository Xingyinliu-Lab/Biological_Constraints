# Public scoring call

```bash
python constraints/protein_sequence/PS04/scoring/score_ps04.py \
  --representations work/PS04/esmc-600m/forward/esmc/embeddings.safetensors \
  --representation-index work/PS04/esmc-600m/forward/esmc/embedding_index.tsv \
  --model-id esmc-600m \
  --output-dir work/PS04/esmc-600m/score
```
