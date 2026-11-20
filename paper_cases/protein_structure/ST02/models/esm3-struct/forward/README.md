# ST02 ESM3-STRUCT Forward

Use the environment and checkpoint identity recorded in [`paper_cases/model_resources/protein_structure/esm3-struct/`](../../../../../model_resources/protein_structure/esm3-struct/README.md). The Forward input is the complete set of 509 released PDB files in `constraints/protein_structure/ST02/refset/coordinates/`.

```bash
python paper_cases/protein_structure/ST02/models/esm3-struct/forward/esm3_structure_mean_extract.py \
  --input-dir constraints/protein_structure/ST02/refset/coordinates \
  --objects constraints/protein_structure/ST02/refset/objects.tsv \
  --output work/ST02/esm3-struct/embeddings.npz
```

The registered structure-only interface masks all interior sequence tokens, encodes structure tokens plus backbone coordinates, and mean-pools non-boundary residue embeddings. The output NPZ contains `object_ids` `(509,)` and `embeddings` `(509,d)`. A `FORWARD_REPORT.json` is written beside it. Score the NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
