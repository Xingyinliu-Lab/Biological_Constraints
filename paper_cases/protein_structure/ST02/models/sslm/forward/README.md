# ST02 SSLM Forward

Use the environment and checkpoint recorded in [`paper_cases/model_resources/protein_structure/sslm/`](../../../../../model_resources/protein_structure/sslm/README.md). SSLM consumes Foldseek 3Di sequences derived from the complete set of 509 released PDB files, not amino-acid sequences.

## Generate the 3Di input

```bash
foldseek createdb constraints/protein_structure/ST02/refset/coordinates work/ST02/sslm/foldseek/DB
foldseek lndb work/ST02/sslm/foldseek/DB_h work/ST02/sslm/foldseek/DB_ss_h
foldseek convert2fasta work/ST02/sslm/foldseek/DB work/ST02/sslm/structures_aa.fasta
foldseek convert2fasta work/ST02/sslm/foldseek/DB_ss work/ST02/sslm/structures_3di.fasta
```

Keep the 3Di FASTA identifiers emitted by Foldseek and corresponding to the
coordinate-file name stems. The encoder uses `--objects` to map those stems to
the formal ST02 `record_id`; its output NPZ contains the formal IDs.

## Encode and pool

```bash
python paper_cases/protein_structure/ST02/models/sslm/forward/sslm_mean_extract.py \
  --fasta work/ST02/sslm/structures_3di.fasta \
  --objects constraints/protein_structure/ST02/refset/objects.tsv \
  --model-dir /path/to/registered/sslm-model \
  --output work/ST02/sslm/embeddings.npz
```

The registered readout is the unweighted mean of non-special 3Di token states. The output NPZ contains `object_ids` `(509,)` and `embeddings` `(509,d)`; `FORWARD_REPORT.json` is written beside it. Score the NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
