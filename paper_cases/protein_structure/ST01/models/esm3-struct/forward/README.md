# ESM3-STRUCT Forward

The extractor directly produces the global representation; it does not use the residue pooling helper.

```bash
mkdir -p work/ST01/esm3-struct
python paper_cases/protein_structure/ST01/models/esm3-struct/forward/scripts/esm3_structure_mean_extract.py --input-dir constraints/protein_structure/ST01/refset/coordinates --objects constraints/protein_structure/ST01/refset/objects.tsv --output work/ST01/esm3-struct/embeddings.npz
```

It uses the structure-only ESM3 input, masks interior sequence tokens, and mean-pools residue states excluding boundaries. `objects.tsv` maps coordinate stems to canonical output IDs. Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md); see [`model_resources/protein_structure/esm3-struct`](../../../../../model_resources/protein_structure/esm3-struct/README.md) for environment identity.
