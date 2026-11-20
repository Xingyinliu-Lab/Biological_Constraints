
# DNA01 × crossdna-8m: execution commands

The registered execution parameters are in `RUNPLAN.json`. The portable Forward entry is:

```bash
python paper_cases/dna/DNA01/models/crossdna-8m/forward/scripts/encode_dna01_crossdna.py \
  --model-dir <CHECKPOINT_DIRECTORY> \
  --object-table constraints/dna/DNA01/refset/DNA01_OBJECT_TABLE.tsv \
  --fasta constraints/dna/DNA01/refset/records.fasta \
  --output work/DNA01/crossdna-8m/forward
```

Replace angle-bracketed model locations with the exact resources identified under `paper_cases/model_resources/dna/crossdna-8m/`. Do not change the RefSet population or extraction semantics.
