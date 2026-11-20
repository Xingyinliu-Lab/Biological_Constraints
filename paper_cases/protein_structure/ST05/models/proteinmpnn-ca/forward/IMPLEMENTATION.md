# proteinmpnn-ca Forward implementation

Use the registered [ProteinMPNN-CA model resource](../../../../../model_resources/protein_structure/proteinmpnn-ca/README.md). From the repository root:

```bash
python paper_cases/protein_structure/ST05/forward_tools/stage_st05_coordinate_panel.py \
  --refset-dir constraints/protein_structure/ST05/refset --output-dir work/ST05/proteinmpnn-ca/coordinates
python paper_cases/protein_structure/ST05/models/proteinmpnn-ca/forward/scripts/encode_proteinmpnn_structure.py \
  --input-dir work/ST05/proteinmpnn-ca/coordinates --output-dir work/ST05/proteinmpnn-ca/residue \
  --proteinmpnn-repo /path/to/ProteinMPNN --weights /path/to/v_48_020.pt --device cuda:0
python paper_cases/protein_structure/ST05/models/proteinmpnn-ca/forward/scripts/pool_st05_residue_embeddings.py \
  --relations constraints/protein_structure/ST05/refset/relations.tsv --residue-dir work/ST05/proteinmpnn-ca/residue \
  --model-id proteinmpnn-ca --output work/ST05/proteinmpnn-ca/forward/embeddings.npz
```

The encoder's CA structural representation is pooled by the registered unweighted global residue mean. Score the NPZ with `../scoring_call/COMMANDS.md`.

Retained scripts:

- `encode_proteinmpnn_structure.py`
- `pool_st05_residue_embeddings.py`
