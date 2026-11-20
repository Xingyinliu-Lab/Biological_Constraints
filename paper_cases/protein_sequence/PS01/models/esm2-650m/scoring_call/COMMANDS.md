# PS01 scoring call for esm2-650m

After regenerating the registered safetensors Forward output, run the unique public scorer:

```bash
python constraints/protein_sequence/PS01/scoring/score_ps01.py \
  --representations work/PS01/esm2-650m/forward/embeddings.safetensors \
  --representation-index work/PS01/esm2-650m/forward/embedding_index.tsv \
  --model-id esm2-650m --output-dir work/PS01/esm2-650m/score
```
