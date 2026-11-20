
# DNA01 × janusdna: execution commands

The registered execution parameters are in `RUNPLAN.json`. The portable Forward entry is:

```bash
python paper_cases/dna/DNA01/models/janusdna/forward/scripts/encode_dna01_janusdna.py \
  --source <JANUSDNA_SOURCE> --config <CONFIG> --checkpoint <CHECKPOINT> \
  --object-table constraints/dna/DNA01/refset/DNA01_OBJECT_TABLE.tsv \
  --fasta constraints/dna/DNA01/refset/records.fasta \
  --output work/DNA01/janusdna/forward
```

Replace angle-bracketed model locations with the exact resources identified under `paper_cases/model_resources/dna/janusdna/`. Do not change the RefSet population or extraction semantics.
