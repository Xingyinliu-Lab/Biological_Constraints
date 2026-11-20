# ST04 esm-if1 scoring call

This case uses the documented ESM model-input limitation; it does not define a second ST04 scoring method.

```bash
python constraints/protein_structure/ST04/scoring/score_st04.py \
  --relations paper_cases/protein_structure/ST04/scope/ESM_INPUT_LIMITED_RELATIONS.tsv \
  --scope-record paper_cases/protein_structure/ST04/scope/ESM_INPUT_LIMITATION.json \
  --refset-release constraints/protein_structure/ST04/refset/ST04_ASSAY_REFSET_RELEASE_v0.1.1_FROZEN_20261012.json \
  --numeric-method constraints/protein_structure/ST04/scoring/ST04_NUMERIC_METHOD_v0.1.4_F_ONLY.json \
  --split-release constraints/protein_structure/ST04/refset/ST04_SPLIT_RELEASE_v0.1.1_FROZEN_20261012.json \
  --embeddings work/ST04/esm-if1/ENDPOINT_EMBEDDINGS.npz \
  --model-id ESM-IF1 \
  --out-dir work/ST04/esm-if1/score
```
