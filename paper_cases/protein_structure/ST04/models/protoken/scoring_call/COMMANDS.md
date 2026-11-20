# ST04 protoken scoring call

The registered paper case uses the frozen FULL relation scope.

```bash
python constraints/protein_structure/ST04/scoring/score_st04.py \
  --relations constraints/protein_structure/ST04/refset/relations.tsv \
  --refset-release constraints/protein_structure/ST04/refset/ST04_ASSAY_REFSET_RELEASE_v0.1.1_FROZEN_20261012.json \
  --numeric-method constraints/protein_structure/ST04/scoring/ST04_NUMERIC_METHOD_v0.1.4_F_ONLY.json \
  --split-release constraints/protein_structure/ST04/refset/ST04_SPLIT_RELEASE_v0.1.1_FROZEN_20261012.json \
  --embeddings work/ST04/protoken/ENDPOINT_EMBEDDINGS.npz \
  --model-id ProToken \
  --out-dir work/ST04/protoken/score
```
