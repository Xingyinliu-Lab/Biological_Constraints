# PS01 scoring call for prott5

After regenerating the registered safetensors Forward output, run the unique public scorer:

```bash
python constraints/protein_sequence/PS01/scoring/score_ps01.py \
  --representations work/PS01/prott5/forward/embeddings.safetensors \
  --representation-index work/PS01/prott5/forward/embedding_index.tsv \
  --model-id prott5 --output-dir work/PS01/prott5/score
```
