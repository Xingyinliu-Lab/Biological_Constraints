# PS01 scoring call for esmc-600m

After regenerating the registered safetensors Forward output, run the unique public scorer:

```bash
python constraints/protein_sequence/PS01/scoring/score_ps01.py \
  --representations work/PS01/esmc-600m/forward/embeddings.safetensors \
  --representation-index work/PS01/esmc-600m/forward/embedding_index.tsv \
  --model-id esmc-600m --output-dir work/PS01/esmc-600m/score
```
