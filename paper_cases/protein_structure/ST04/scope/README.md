# ST04 paper-case scope assets

`ESM_INPUT_LIMITED_RELATIONS.tsv` is the documented scoring input for ESM-IF1 and ESM3-STRUCT. It contains 152,173 relations, 73,565 directed queries and 48,279 residue endpoints. `ESM_INPUT_LIMITATION.json` records the two structural-input exclusions, the three resulting query exclusions and all seven excluded relation IDs. The remaining three paper models use the complete public relation table at `constraints/protein_structure/ST04/refset/relations.tsv`.

`score_st04.py` scores the complete public RefSet by default. A non-FULL input requires `--scope-record`, whose relation-table hash and population counts are checked before scoring. This preserves one RefSet and one scoring method while documenting an unavoidable model-input limitation.

`FORMAL_ENDPOINT_TOKEN_MAP.tsv` maps each formal residue endpoint to the token index produced by the registered model implementation. Its columns are:

- `model_id`: registered model identity;
- `residue_id`: formal ST04 endpoint identity;
- `structure_object_id`: public RefSet structure identity;
- `token_index`: zero-based index in that model's residue-level Forward output.

The map is shared execution evidence, not a second RefSet or a scoring algorithm. It contains only the endpoint rows required by the complete RefSet or, for the two ESM cases, their documented model-input limitation. `FORMAL_ENDPOINT_TOKEN_MAP_IDENTITY.json` records the released map hash and endpoint counts. The endpoint materializer checks the released map against `constraints/protein_structure/ST04/refset/residue_evidence.tsv` before producing scorer input.
