# ProToken Forward

Use upstream ProToken stage-1 inference with its encoder, decoder and VQ configuration files.

```bash
cd /path/to/ProToken/PROTOKEN
PYTHONPATH=. python scripts/infer_batch.py --encoder_config config/encoder.yaml --decoder_config config/decoder.yaml --vq_config config/vq.yaml --pdb_dir_path /repo/constraints/protein_structure/ST01/refset/coordinates --save_dir_path /repo/work/ST01/protoken/stage1 --load_ckpt_path /path/to/protoken_params_100000.pkl
cd /repo
python paper_cases/protein_structure/ST01/models/protoken/forward/scripts/export_protoken_stage1.py --token-index work/ST01/protoken/stage1/stage_1/vq_code_indexes_dict.pkl --checkpoint /path/to/protoken_params_100000.pkl --objects constraints/protein_structure/ST01/refset/objects.tsv --output work/ST01/protoken/embeddings.npz
```

The exporter maps stage-1 token IDs through the checkpoint VQ codebook and writes canonical ST01 IDs. Then run [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
