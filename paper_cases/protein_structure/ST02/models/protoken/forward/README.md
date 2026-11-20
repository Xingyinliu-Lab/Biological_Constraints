# ST02 ProToken Forward

Use the ProToken source, configuration, environment and checkpoint recorded in [`paper_cases/model_resources/protein_structure/protoken/`](../../../../../model_resources/protein_structure/protoken/README.md). The Forward input is the complete set of 509 released PDB files in `constraints/protein_structure/ST02/refset/coordinates/`.

## Generate structure tokens

Stage copied inputs first. This leaves the released RefSet untouched and only
renumbers residue labels sequentially within a chain for structures that carry
PDB insertion codes, which the registered upstream ProToken source rejects.

```bash
python paper_cases/protein_structure/ST02/models/protoken/forward/prepare_protoken_inputs.py \
  --objects constraints/protein_structure/ST02/refset/objects.tsv \
  --input-dir constraints/protein_structure/ST02/refset/coordinates \
  --output-dir work/ST02/protoken/coordinates
```

Run the frozen ProToken stage-1 inference from the registered upstream source checkout:

```bash
python paper_cases/protein_structure/ST02/models/protoken/forward/infer_batch_st02_batch1.py \
  --encoder_config /path/to/ProToken/config/encoder.yaml \
  --decoder_config /path/to/ProToken/config/decoder.yaml \
  --vq_config /path/to/ProToken/config/vq.yaml \
  --pdb_dir_path work/ST02/protoken/coordinates \
  --save_dir_path work/ST02/protoken/raw \
  --load_ckpt_path /path/to/protoken_params_100000.pkl \
  --padding_len 1300 \
  --random_seed 8888 \
  --np_random_seed 18888
```

The required token dictionary is `work/ST02/protoken/raw/stage_1/vq_code_indexes_dict.pkl`. Convert the registered code IDs through the same checkpoint codebook and apply the frozen unweighted token mean:

```bash
python paper_cases/protein_structure/ST02/models/protoken/forward/finish_protoken_st02.py \
  --tokens work/ST02/protoken/raw/stage_1/vq_code_indexes_dict.pkl \
  --checkpoint /path/to/protoken_params_100000.pkl \
  --objects constraints/protein_structure/ST02/refset/objects.tsv \
  --relations constraints/protein_structure/ST02/refset/relations.tsv \
  --output work/ST02/protoken/embeddings.npz
```

The final NPZ contains `object_ids` `(509,)` and `embeddings` `(509,d)`. Score it with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
