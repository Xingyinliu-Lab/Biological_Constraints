# sslm Forward implementation

Use the registered [SSLM model resource](../../../../../model_resources/protein_structure/sslm/README.md). From the repository root, stage the released structures, then derive the required Foldseek 3Di FASTA with record IDs unchanged:

```bash
python paper_cases/protein_structure/ST05/forward_tools/stage_st05_coordinate_panel.py \
  --refset-dir constraints/protein_structure/ST05/refset --output-dir work/ST05/sslm/coordinates
foldseek createdb work/ST05/sslm/coordinates work/ST05/sslm/foldseek_db
foldseek convert2fasta work/ST05/sslm/foldseek_db_ss work/ST05/sslm/structures_3di.fasta
python paper_cases/protein_structure/ST05/models/sslm/forward/scripts/sslm_residue_extract_st05.py \
  --fasta work/ST05/sslm/structures_3di.fasta --model-dir /path/to/sslm-model --output-dir work/ST05/sslm/residue
python paper_cases/protein_structure/ST05/models/sslm/forward/scripts/pool_st05_residue_embeddings.py \
  --relations constraints/protein_structure/ST05/refset/relations.tsv --residue-dir work/ST05/sslm/residue \
  --model-id sslm --output work/ST05/sslm/forward/embeddings.npz
```

The readout is the retained global mean of non-special 3Di-token states. Score the NPZ with `../scoring_call/COMMANDS.md`.

Retained scripts:

- `pool_st05_residue_embeddings.py`
- `sslm_residue_extract_st05.py`
