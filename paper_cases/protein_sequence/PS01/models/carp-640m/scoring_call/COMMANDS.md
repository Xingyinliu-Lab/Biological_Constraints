# PS01 scoring call for carp-640m

After regenerating the registered safetensors Forward output, run the unique public scorer:

```bash
python constraints/protein_sequence/PS01/scoring/score_ps01.py \
  --representations work/PS01/carp-640m/forward/embeddings.safetensors \
  --representation-index work/PS01/carp-640m/forward/embedding_index.tsv \
  --model-id carp-640m --output-dir work/PS01/carp-640m/score
```
