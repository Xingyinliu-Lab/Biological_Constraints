# PS09 public scoring call

```bash
python constraints/protein_sequence/PS09/scoring/score_ps09.py \
  --representations work/PS09/esm2-650m/forward/endpoint_embeddings.safetensors \
  --model-id esm2-650m \
  --output-dir work/PS09/esm2-650m/score
```

This command uses the released PS09 RefSet, frozen five-fold readout and public aggregation. Copy the completed output directory into this case's `output/` when registering the result.
