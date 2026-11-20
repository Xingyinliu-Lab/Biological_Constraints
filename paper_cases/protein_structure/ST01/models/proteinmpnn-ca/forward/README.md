# ProteinMPNN-CA Forward

From the repository root:

```bash
mkdir -p work/ST01/proteinmpnn-ca/residue
python paper_cases/protein_structure/ST01/models/proteinmpnn-ca/forward/scripts/encode_proteinmpnn_structure.py --input-dir constraints/protein_structure/ST01/refset/coordinates --output-dir work/ST01/proteinmpnn-ca/residue --proteinmpnn-repo /path/to/ProteinMPNN --weights /path/to/v_48_020.pt
python paper_cases/protein_structure/ST01/shared_forward/pool_residue_embeddings.py --objects constraints/protein_structure/ST01/refset/objects.tsv --residue-dir work/ST01/proteinmpnn-ca/residue --output work/ST01/proteinmpnn-ca/embeddings.npz
```

The shared pooling step maps coordinate filenames to canonical ST01 object IDs. Score the resulting NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md). See [`model_resources/protein_structure/proteinmpnn-ca`](../../../../../model_resources/protein_structure/proteinmpnn-ca/README.md) for the environment and weights.
