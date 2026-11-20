# Public scoring call

```bash
python constraints/protein_sequence/PS04/scoring/score_ps04.py \
  --representations work/PS04/prott5/forward/prott5/embeddings.safetensors \
  --representation-index work/PS04/prott5/forward/prott5/embedding_index.tsv \
  --model-id prott5 \
  --output-dir work/PS04/prott5/score
```
