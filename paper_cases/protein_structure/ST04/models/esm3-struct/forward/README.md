# ST04 ESM3-STRUCT Forward

Use the ESM3-STRUCT environment and checkpoint recorded in `paper_cases/model_resources/protein_structure/esm3-struct/`. This paper case has a documented ESM structural-input limitation; see `scope/ESM_INPUT_LIMITATION.json`.

## Inputs

- the 586 formal PDB files in `constraints/protein_structure/ST04/refset/objects/experimental_pdb/`;
- the ESM input-limited relation table, endpoint token map and residue evidence listed in `paper_cases/protein_structure/ST04/README.md`;
- the registered ESM3-STRUCT checkpoint.

## Forward and endpoint materialization

```bash
python paper_cases/protein_structure/ST04/models/esm3-struct/forward/esm3_structure_residue_extract_st04_v0_1_2.py \
  --input-dir constraints/protein_structure/ST04/refset/objects/experimental_pdb \
  --output-dir work/ST04/esm3-struct/residue \
  --checkpoint-binding 5ead5a135c658068db6a4f1b933e72d6110992c4668822e1c0e2dcc53e38acd9

python paper_cases/protein_structure/ST04/models/esm3-struct/forward/materialize_st04_endpoint_embeddings_v0_1_2.py \
  --model-id ESM3-STRUCT \
  --relations paper_cases/protein_structure/ST04/scope/ESM_INPUT_LIMITED_RELATIONS.tsv \
  --token-map paper_cases/protein_structure/ST04/scope/FORMAL_ENDPOINT_TOKEN_MAP.tsv \
  --residue-evidence constraints/protein_structure/ST04/refset/residue_evidence.tsv \
  --residue-dir work/ST04/esm3-struct/residue \
  --output work/ST04/esm3-struct/ENDPOINT_EMBEDDINGS.npz
```

The materialized NPZ contains `residue_ids` `(n,)` and `embeddings` `(n,d)`. Score it with `../scoring_call/COMMANDS.md`.
