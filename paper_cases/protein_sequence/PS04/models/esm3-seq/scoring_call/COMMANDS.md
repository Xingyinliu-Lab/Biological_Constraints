# Public scoring call

```bash
python constraints/protein_sequence/PS04/scoring/score_ps04.py \
  --representations work/PS04/esm3-seq/forward/esm3/embeddings.safetensors \
  --representation-index work/PS04/esm3-seq/forward/esm3/embedding_index.tsv \
  --model-id esm3-seq \
  --output-dir work/PS04/esm3-seq/score
```
