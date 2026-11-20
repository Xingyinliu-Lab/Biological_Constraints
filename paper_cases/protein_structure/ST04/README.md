# ST04 paper implementation cases

These five cases preserve the paper Forward implementations and registered F-only outputs. They do not contain another scoring implementation: every case hands `ENDPOINT_EMBEDDINGS.npz` to the public ST04 scorer in `constraints/protein_structure/ST04/scoring/`.

## Shared assets

- Formal coordinate input: `constraints/protein_structure/ST04/refset/objects/experimental_pdb/` (586 admitted structures).
- Frozen residue mapping: `constraints/protein_structure/ST04/refset/residue_evidence.tsv`.
- Formal endpoint-to-model-token mapping: `scope/FORMAL_ENDPOINT_TOKEN_MAP.tsv`.
- Complete formal RefSet: `constraints/protein_structure/ST04/refset/relations.tsv`.
- ESM input-limited relation table: `scope/ESM_INPUT_LIMITED_RELATIONS.tsv`, documented by `scope/ESM_INPUT_LIMITATION.json`.

The released token map is the exact endpoint closure of the historical executed map: 241,404 model-endpoint rows, SHA256 `0d2402cffa6ba83f7a956f128b51de24a4857c8ce1c03f23c12483a5dd1c5e5c`. Its derivation and source identity are recorded in `scope/FORMAL_ENDPOINT_TOKEN_MAP_IDENTITY.json`. The residue evidence has 98,209 data rows and SHA256 `e0d1c4fa106959fa4ee1c142f0ec18576b806fcb8c7cf45f4a35df300e653c05`.

## Cases

| Model | Scoring input | Forward manual |
|---|---|---|
| ESM-IF1 | 152,173 relations; documented ESM input limitation | `models/esm-if1/forward/README.md` |
| ESM3-STRUCT | 152,173 relations; documented ESM input limitation | `models/esm3-struct/forward/README.md` |
| ProteinMPNN-CA | complete formal RefSet (152,180 relations) | `models/proteinmpnn-ca/forward/README.md` |
| ProToken | complete formal RefSet (152,180 relations) | `models/protoken/forward/README.md` |
| SSLM | complete formal RefSet (152,180 relations) | `models/sslm/forward/README.md` |

Each manual identifies the exact extractor, its inputs and outputs, the endpoint materialization step and the public scoring handoff. Model environments and checkpoint identities are maintained once under `paper_cases/model_resources/protein_structure/`.

Formal embedding NPZ files, fit arrays, OOF predictions and relation Ledgers are reproducible process outputs and are not distributed in the public Paper Cases package. Public reproduction begins with the released RefSet and model Forward manual, then runs the public scorer to regenerate the complete F audit chain. The published `FORMAL_RESULT.json` files are score summaries, not claims of byte-level identity to non-public historical intermediates.
