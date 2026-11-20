# ST05

Use this page to evaluate a model against the ST05 constraint. The commands below connect the released RefSet, representation interface and public F-only scorer.

## Public scoring interface

- Official public entry: `constraints/protein_structure/ST05/scoring/score_st05.py`
- Scoring environment: [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md)
- All compatible representations use the same frozen RefSet, numeric method, aggregation and uncertainty procedure.
## Run the scorer

### Run the official invariant scorer

Provide an NPZ containing exactly `base_object_ids`, `relation_ids`,
`base_embeddings`, and `variant_embeddings`. The public entry uses the
frozen invariant law and derives the implementation-precision tolerance from
the declared output dtype and representation dimension:

`base_object_ids` and `relation_ids` are unique one-dimensional string arrays.
`base_embeddings` has shape `(number_of_base_objects, d)` and
`variant_embeddings` has shape `(number_of_relations, d)` with the same
positive `d`; rows align to their corresponding ID arrays. Both ID sets must
exactly equal the released relation universe, and all values must be finite.

Both stored embedding arrays must have exactly the dtype supplied as
`--output-dtype`; the scorer verifies this before its internal float64
calculation. The declared dtype, rather than internal Forward autocast dtype,
sets the prespecified tolerance precision.

```bash
python constraints/protein_structure/ST05/scoring/score_st05.py \
  --relations constraints/protein_structure/ST05/refset/relations.tsv \
  --refset-release constraints/protein_structure/ST05/refset/ST05_ASSAY_REFSET_RELEASE_v0.1.0_FROZEN_20261012.json \
  --embeddings my_model_st05_embeddings.npz \
  --model-id my_model \
  --output-dtype float32 \
  --out work/ST05/my_model/EVALUATION_RESULT.json
```

This entry implements the invariant branch used by all five published paper
models. A future third-party model claiming the separately defined executable
equivariant branch must preregister and validate its action operators; it must
not use this invariant entry.

`EVALUATION_RESULT.json` reports `F`, `F_95_CI`, invariant-branch identity,
population counts, derived numerical tolerances and residual summaries. ST05
uses the invariant satisfaction fraction directly as F; it does not introduce
a chance-centered transform. The scorer
also writes the per-transform ledger used to derive these summaries.

## Formal method files

Frozen scoring implementation:

- `constraints/protein_structure/ST05/scoring/st05_f_only_core.py`
- `constraints/protein_structure/ST05/scoring/ST05_NUMERIC_METHOD_v0.1.2_F_ONLY.json`

## Released RefSet

Ready-to-use RefSet assets:

- `constraints/protein_structure/ST05/refset/attrition.tsv`
- `constraints/protein_structure/ST05/refset/dependency.tsv`
- `constraints/protein_structure/ST05/refset/manifest.json`
- `constraints/protein_structure/ST05/refset/objects.tsv`
- `constraints/protein_structure/ST05/refset/relations.tsv`
- `constraints/protein_structure/ST05/refset/ST05_ASSAY_REFSET_RELEASE_v0.1.0_FROZEN_20261012.json`
- `constraints/protein_structure/ST05/refset/transforms.tsv`
