# ST02 ProteinMPNN-CA Forward

Use the environment, source revision and checkpoint recorded in [`paper_cases/model_resources/protein_structure/proteinmpnn-ca/`](../../../../../model_resources/protein_structure/proteinmpnn-ca/README.md). The Forward input is the complete set of 509 released PDB files in `constraints/protein_structure/ST02/refset/coordinates/`.

## Generate residue representations

```bash
python paper_cases/protein_structure/ST02/models/proteinmpnn-ca/forward/encode_proteinmpnn_structure.py \
  --input-dir constraints/protein_structure/ST02/refset/coordinates \
  --output-dir work/ST02/proteinmpnn-ca/residue \
  --proteinmpnn-repo /path/to/registered/ProteinMPNN \
  --weights /path/to/registered/v_48_020.pt \
  --device cuda:0
```

The registered path uses the CA-only encoder and retains `residue_embedding` `(L,128)` for every object. Apply the frozen unweighted residue mean:

```bash
python paper_cases/protein_structure/ST02/forward_tools/pool_residue_embeddings.py \
  --objects constraints/protein_structure/ST02/refset/objects.tsv \
  --relations constraints/protein_structure/ST02/refset/relations.tsv \
  --residue-dir work/ST02/proteinmpnn-ca/residue \
  --output work/ST02/proteinmpnn-ca/embeddings.npz
```

The final NPZ contains `object_ids` `(509,)` and `embeddings` `(509,128)`. Score it with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
