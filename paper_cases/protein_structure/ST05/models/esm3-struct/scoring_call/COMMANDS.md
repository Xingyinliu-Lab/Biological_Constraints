# ST05 scoring call for esm3-struct

The Forward implementation must produce the registered four-array NPZ. Run the unique public scorer:

```bash
python constraints/protein_structure/ST05/scoring/score_st05.py \
  --relations constraints/protein_structure/ST05/refset/relations.tsv \
  --refset-release constraints/protein_structure/ST05/refset/ST05_ASSAY_REFSET_RELEASE_v0.1.0_FROZEN_20261012.json \
  --embeddings work/ST05/esm3-struct/forward/embeddings.npz \
  --model-id esm3-struct --output-dtype float32 \
  --out work/ST05/esm3-struct/score/EVALUATION_RESULT.json
```

The scorer also writes `ST05_TRANSFORM_RESULT_LEDGER.tsv` beside the JSON result.
