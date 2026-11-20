# ST04 ESM-IF1 Forward

Use the ESM-IF1 environment and checkpoint recorded in `paper_cases/model_resources/protein_structure/esm-if1/`. This paper case has a documented ESM structural-input limitation; see `scope/ESM_INPUT_LIMITATION.json`.

## Inputs

- `constraints/protein_structure/ST04/refset/objects/experimental_pdb/`: 586 formal PDB inputs;
- `paper_cases/protein_structure/ST04/scope/ESM_INPUT_LIMITED_RELATIONS.tsv`;
- `paper_cases/protein_structure/ST04/scope/FORMAL_ENDPOINT_TOKEN_MAP.tsv`;
- `constraints/protein_structure/ST04/refset/residue_evidence.tsv`;
- an installed ESM repository and the registered ESM-IF1 checkpoint.

## Forward

```bash
python paper_cases/protein_structure/ST04/models/esm-if1/forward/encode_esm_if1_incremental_v0_1_2.py \
  --input-dir constraints/protein_structure/ST04/refset/objects/experimental_pdb \
  --output-dir work/ST04/esm-if1/residue \
  --esm-repo /path/to/esm \
  --checkpoint-binding be4ba36edec22a9bfaa4946ff6b2815f1f19d8a3d7e0eada8b796d5a0eae9fd4
```

The extractor writes one NPZ per structure containing residue-level representations plus `incremental_encoding_summary.json`.

## Materialize scorer input

```bash
python paper_cases/protein_structure/ST04/models/esm-if1/forward/materialize_st04_endpoint_embeddings_v0_1_2.py \
  --model-id ESM-IF1 \
  --relations paper_cases/protein_structure/ST04/scope/ESM_INPUT_LIMITED_RELATIONS.tsv \
  --token-map paper_cases/protein_structure/ST04/scope/FORMAL_ENDPOINT_TOKEN_MAP.tsv \
  --residue-evidence constraints/protein_structure/ST04/refset/residue_evidence.tsv \
  --residue-dir work/ST04/esm-if1/residue \
  --output work/ST04/esm-if1/ENDPOINT_EMBEDDINGS.npz
```

The final NPZ contains `residue_ids` `(n,)` and `embeddings` `(n,d)`. It is the input consumed by the command in `../scoring_call/COMMANDS.md`.
