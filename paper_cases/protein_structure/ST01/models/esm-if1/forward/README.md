# ESM-IF1 Forward

From the repository root, encode the 28 released PDBs, then pool each residue array by the canonical ID in `objects.tsv`:

```bash
mkdir -p work/ST01/esm-if1/residue
python paper_cases/protein_structure/ST01/models/esm-if1/forward/scripts/encode_esm_if1_incremental.py --input-dir constraints/protein_structure/ST01/refset/coordinates --output-dir work/ST01/esm-if1/residue --esm-repo /path/to/esm
python paper_cases/protein_structure/ST01/shared_forward/pool_residue_embeddings.py --objects constraints/protein_structure/ST01/refset/objects.tsv --residue-dir work/ST01/esm-if1/residue --output work/ST01/esm-if1/embeddings.npz
```

`embeddings.npz` contains canonical `object_ids` and whole-structure residue-mean `embeddings`. Run [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md) next. The required checkpoint and environment are in [`model_resources/protein_structure/esm-if1`](../../../../../model_resources/protein_structure/esm-if1/README.md).
