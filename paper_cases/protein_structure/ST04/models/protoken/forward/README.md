# ST04 ProToken Forward

Use the ProToken environment, source and checkpoint recorded in `paper_cases/model_resources/protein_structure/protoken/`. This case uses the complete FULL relation population.

## Prepare structures and run ProToken

The registered implementation normalizes PDB insertion-code fields without changing atom identity, then runs the frozen ProToken encoder configuration.

Install the registered ProToken source tree, then expose its repository root before running the encoder; `infer_batch_st02_batch1.py` imports its `model` and `tokenizer` packages directly:

```powershell
$env:PYTHONPATH = "<PROTOKEN_SOURCE_ROOT>$([IO.Path]::PathSeparator)$env:PYTHONPATH"
```

```bash
python paper_cases/protein_structure/ST04/models/protoken/forward/prepare_protoken_inputs_st04.py \
  --input-dir constraints/protein_structure/ST04/refset/objects/experimental_pdb \
  --output-dir work/ST04/protoken/input_pdb

python paper_cases/protein_structure/ST04/models/protoken/forward/infer_batch_st02_batch1.py \
  --encoder_config /path/to/ProToken/config/encoder.yaml \
  --decoder_config /path/to/ProToken/config/decoder.yaml \
  --vq_config /path/to/ProToken/config/vq.yaml \
  --pdb_dir_path work/ST04/protoken/input_pdb \
  --save_dir_path work/ST04/protoken/raw \
  --load_ckpt_path /path/to/protoken_params_100000.pkl \
  --padding_len 1300 \
  --random_seed 8888 \
  --np_random_seed 18888
```

The required token dictionary is `work/ST04/protoken/raw/stage_1/vq_code_indexes_dict.pkl`.

## Materialize scorer input

```bash
python paper_cases/protein_structure/ST04/models/protoken/forward/materialize_st04_endpoint_embeddings_v0_1_2.py \
  --model-id ProToken \
  --relations constraints/protein_structure/ST04/refset/relations.tsv \
  --token-map paper_cases/protein_structure/ST04/scope/FORMAL_ENDPOINT_TOKEN_MAP.tsv \
  --residue-evidence constraints/protein_structure/ST04/refset/residue_evidence.tsv \
  --protoken-tokens work/ST04/protoken/raw/stage_1/vq_code_indexes_dict.pkl \
  --protoken-checkpoint /path/to/protoken_params_100000.pkl \
  --output work/ST04/protoken/ENDPOINT_EMBEDDINGS.npz
```

The output contains `residue_ids` `(n,)` and `embeddings` `(n,d)` and is scored with `../scoring_call/COMMANDS.md`.
