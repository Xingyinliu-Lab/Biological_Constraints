# ST02 ESM-IF1 Forward

Use the environment, source and checkpoint identity recorded in [`paper_cases/model_resources/protein_structure/esm-if1/`](../../../../../model_resources/protein_structure/esm-if1/README.md). The Forward input is the complete set of 509 released PDB files in `constraints/protein_structure/ST02/refset/coordinates/`.

## Generate residue representations

```bash
python paper_cases/protein_structure/ST02/models/esm-if1/forward/encode_esm_if1_incremental.py \
  --input-dir constraints/protein_structure/ST02/refset/coordinates \
  --output-dir work/ST02/esm-if1/residue \
  --esm-repo /path/to/registered/esm-repository
```

Each object file contains `residue_embedding` with shape `(L,d)`. The registered ST02 readout is the unweighted arithmetic mean over its `L` residue rows:

```bash
python paper_cases/protein_structure/ST02/forward_tools/pool_residue_embeddings.py \
  --objects constraints/protein_structure/ST02/refset/objects.tsv \
  --relations constraints/protein_structure/ST02/refset/relations.tsv \
  --residue-dir work/ST02/esm-if1/residue \
  --output work/ST02/esm-if1/embeddings.npz
```

The final NPZ contains canonical ST02 `object_ids` `(509,)` and `embeddings` `(509,d)`. Score it with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md). The scorer, not this Forward code, computes the frozen contrast relations, hierarchy, ID-level scores and aggregate statistics.
