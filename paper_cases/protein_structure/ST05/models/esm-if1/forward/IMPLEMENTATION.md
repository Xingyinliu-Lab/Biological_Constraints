# esm-if1 Forward implementation

Use the registered [ESM-IF1 model resource](../../../../../model_resources/protein_structure/esm-if1/README.md). From the repository root, flatten the released panel without changing coordinate bytes:

```bash
python paper_cases/protein_structure/ST05/forward_tools/stage_st05_coordinate_panel.py \
  --refset-dir constraints/protein_structure/ST05/refset \
  --output-dir work/ST05/esm-if1/coordinates
python paper_cases/protein_structure/ST05/models/esm-if1/forward/scripts/encode_esm_if1_incremental.py \
  --input-dir work/ST05/esm-if1/coordinates --output-dir work/ST05/esm-if1/residue \
  --esm-repo /path/to/esm
python paper_cases/protein_structure/ST05/models/esm-if1/forward/scripts/pool_st05_residue_embeddings.py \
  --relations constraints/protein_structure/ST05/refset/relations.tsv --residue-dir work/ST05/esm-if1/residue \
  --model-id esm-if1 --output work/ST05/esm-if1/forward/embeddings.npz
```

The encoder retains N/CA/C complete residues; pooling is an unweighted global residue mean (float64 accumulation, float32 output). Score the resulting NPZ with `../scoring_call/COMMANDS.md`.

Retained scripts:

- `encode_esm_if1_incremental.py`
- `pool_st05_residue_embeddings.py`
