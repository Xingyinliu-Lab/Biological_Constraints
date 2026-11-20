# Public scoring call

```bash
python constraints/protein_sequence/PS04/scoring/score_ps04.py \
  --representations work/PS04/carp-640m/forward/carp/embeddings.safetensors \
  --representation-index work/PS04/carp-640m/forward/carp/embedding_index.tsv \
  --model-id carp-640m \
  --output-dir work/PS04/carp-640m/score
```
