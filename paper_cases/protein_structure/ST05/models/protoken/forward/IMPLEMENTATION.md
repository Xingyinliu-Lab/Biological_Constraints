# protoken Forward implementation

Use the registered [ProToken model resource](../../../../../model_resources/protein_structure/protoken/README.md) and its upstream source checkout. The stage-1 script imports upstream `model/` and `tokenizer/` packages; either install that checkout or export `PYTHONPATH=/path/to/ProToken` before running. From the repository root:

```bash
python paper_cases/protein_structure/ST05/forward_tools/stage_st05_coordinate_panel.py \
  --refset-dir constraints/protein_structure/ST05/refset --output-dir work/ST05/protoken/coordinates
python paper_cases/protein_structure/ST05/models/protoken/forward/scripts/prepare_protoken_inputs_st05.py \
  --input-dir work/ST05/protoken/coordinates --output-dir work/ST05/protoken/normalized
PYTHONPATH=/path/to/ProToken python paper_cases/protein_structure/ST05/models/protoken/forward/scripts/infer_batch_st02_batch1.py \
  --encoder_config /path/to/ProToken/config/encoder.yaml --decoder_config /path/to/ProToken/config/decoder.yaml \
  --vq_config /path/to/ProToken/config/vq.yaml --pdb_dir_path work/ST05/protoken/normalized \
  --save_dir_path work/ST05/protoken/raw --load_ckpt_path /path/to/protoken_params_100000.pkl
python paper_cases/protein_structure/ST05/models/protoken/forward/scripts/finish_protoken_st05.py \
  --relations constraints/protein_structure/ST05/refset/relations.tsv \
  --tokens work/ST05/protoken/raw/stage_1/vq_code_indexes_dict.pkl \
  --checkpoint /path/to/protoken_params_100000.pkl --output work/ST05/protoken/forward/embeddings.npz
```

The final readout is the retained float64 mean of codebook rows, written as float32. Score the NPZ with `../scoring_call/COMMANDS.md`.

Retained scripts:

- `finish_protoken_st05.py`
- `infer_batch_st02_batch1.py`
- `prepare_protoken_inputs_st05.py`
