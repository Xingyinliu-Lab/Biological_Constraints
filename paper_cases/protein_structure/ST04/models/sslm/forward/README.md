# ST04 SSLM Forward

Use the SSLM environment and checkpoint recorded in `paper_cases/model_resources/protein_structure/sslm/`. This case uses the complete FULL relation population.

## Generate the formal 3Di input

SSLM consumes Foldseek 3Di sequences, not amino-acid sequences. Generate them deterministically from the 586 formal PDB files:

```bash
mkdir -p work/ST04/sslm/foldseek
foldseek createdb constraints/protein_structure/ST04/refset/objects/experimental_pdb work/ST04/sslm/foldseek/DB
foldseek lndb work/ST04/sslm/foldseek/DB_h work/ST04/sslm/foldseek/DB_ss_h
foldseek convert2fasta work/ST04/sslm/foldseek/DB work/ST04/sslm/structures_aa.fasta
foldseek convert2fasta work/ST04/sslm/foldseek/DB_ss work/ST04/sslm/structures_3di.fasta
```

The `structures_3di.fasta` record IDs must remain identical to the public `structure_object_id` values.

## Forward and endpoint materialization

```bash
python paper_cases/protein_structure/ST04/models/sslm/forward/sslm_residue_extract_st04_v0_1_2.py \
  --fasta work/ST04/sslm/structures_3di.fasta \
  --model-dir /path/to/registered/sslm-model \
  --output-dir work/ST04/sslm/residue

python paper_cases/protein_structure/ST04/models/sslm/forward/materialize_st04_endpoint_embeddings_v0_1_2.py \
  --model-id SSLM \
  --relations constraints/protein_structure/ST04/refset/relations.tsv \
  --token-map paper_cases/protein_structure/ST04/scope/FORMAL_ENDPOINT_TOKEN_MAP.tsv \
  --residue-evidence constraints/protein_structure/ST04/refset/residue_evidence.tsv \
  --residue-dir work/ST04/sslm/residue \
  --output work/ST04/sslm/ENDPOINT_EMBEDDINGS.npz
```

The materialized file contains `residue_ids` `(n,)` and `embeddings` `(n,d)`. Score it with `../scoring_call/COMMANDS.md`.
