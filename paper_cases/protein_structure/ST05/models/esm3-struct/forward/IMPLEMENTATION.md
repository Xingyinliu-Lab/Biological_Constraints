# esm3-struct Forward implementation

Use the registered [ESM3-structure model resource](../../../../../model_resources/protein_structure/esm3-struct/README.md). From the repository root:

```bash
python paper_cases/protein_structure/ST05/forward_tools/stage_st05_coordinate_panel.py \
  --refset-dir constraints/protein_structure/ST05/refset --output-dir work/ST05/esm3-struct/coordinates
python paper_cases/protein_structure/ST05/models/esm3-struct/forward/scripts/esm3_structure_residue_extract_st05.py \
  --input-dir work/ST05/esm3-struct/coordinates --output-dir work/ST05/esm3-struct/residue
python paper_cases/protein_structure/ST05/models/esm3-struct/forward/scripts/pool_st05_residue_embeddings.py \
  --relations constraints/protein_structure/ST05/refset/relations.tsv --residue-dir work/ST05/esm3-struct/residue \
  --model-id esm3-struct --output work/ST05/esm3-struct/forward/embeddings.npz
```

The structure-only extractor's residue states are globally mean-pooled using the retained script. Score the NPZ with `../scoring_call/COMMANDS.md`.

Retained scripts:

- `esm3_structure_residue_extract_st05.py`
- `pool_st05_residue_embeddings.py`
