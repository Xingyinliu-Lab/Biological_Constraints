# PS01 scoring call for esm3-seq

After regenerating the registered safetensors Forward output, run the unique public scorer:

```bash
python constraints/protein_sequence/PS01/scoring/score_ps01.py \
  --representations work/PS01/esm3-seq/forward/embeddings.safetensors \
  --representation-index work/PS01/esm3-seq/forward/embedding_index.tsv \
  --model-id esm3-seq --output-dir work/PS01/esm3-seq/score
```
