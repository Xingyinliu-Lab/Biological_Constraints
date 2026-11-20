# RNA01-02 scoring

`score_rna0102_global_ci.py` is the only formal public scorer. Both paper-model and third-party representations use it.

## Inputs

- one complete L16K representation NPZ;
- two repeat NPZs for the frozen 128-ladder calibration panel;
- a model ID and output directory.

All files use the exact `sequence_sha256` + `vectors` schema in [`../README.md`](../README.md). The scorer rejects missing, duplicate, extra, complex, Boolean, non-numeric, non-finite or zero-dimensional input. It validates the released sequence objects, edit relations, population membership and clusters from their biological content before scoring.

## Fixed computation

The scorer derives contrast tolerance from repeat calibration, then computes cosine distances, left/right tie-aware credits, ladder scores, outer-cluster means and a 20,000-draw outer-cluster bootstrap. These formal choices are not CLI options.

```bash
python constraints/rna/RNA01-02/scoring/score_rna0102_global_ci.py \
  --vectors <L16K_VECTORS.npz> \
  --calibration-repeat-1 <CALIBRATION_REPEAT_1.npz> \
  --calibration-repeat-2 <CALIBRATION_REPEAT_2.npz> \
  --model-id <MODEL_ID> \
  --output <OUTPUT_DIRECTORY>
```

See [`ENVIRONMENT.md`](ENVIRONMENT.md) for dependencies and [`../README.md`](../README.md) for output columns and interpretation.
