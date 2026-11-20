# ST04 ProteinMPNN-CA Forward

Use the environment, source revision and CA-only checkpoint recorded in `paper_cases/model_resources/protein_structure/proteinmpnn-ca/`. This case uses the complete FULL relation population.

## Forward and endpoint materialization

```bash
python paper_cases/protein_structure/ST04/models/proteinmpnn-ca/forward/encode_proteinmpnn_structure_v0_1_2.py \
  --input-dir constraints/protein_structure/ST04/refset/objects/experimental_pdb \
  --output-dir work/ST04/proteinmpnn-ca/residue \
  --proteinmpnn-repo /path/to/ProteinMPNN \
  --weights /path/to/v_48_020.pt \
  --device cuda:0

python paper_cases/protein_structure/ST04/models/proteinmpnn-ca/forward/materialize_st04_endpoint_embeddings_v0_1_2.py \
  --model-id ProteinMPNN-CA \
  --relations constraints/protein_structure/ST04/refset/relations.tsv \
  --token-map paper_cases/protein_structure/ST04/scope/FORMAL_ENDPOINT_TOKEN_MAP.tsv \
  --residue-evidence constraints/protein_structure/ST04/refset/residue_evidence.tsv \
  --residue-dir work/ST04/proteinmpnn-ca/residue \
  --output work/ST04/proteinmpnn-ca/ENDPOINT_EMBEDDINGS.npz
```

The extractor writes residue-level NPZ files. The materializer validates the frozen ID mapping and writes `residue_ids` `(n,)` plus `embeddings` `(n,d)`. Score that file with `../scoring_call/COMMANDS.md`.
