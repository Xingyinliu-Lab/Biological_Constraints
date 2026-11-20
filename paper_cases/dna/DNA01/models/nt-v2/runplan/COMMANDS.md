
# DNA01 × nt-v2: execution commands

The registered execution parameters are in `RUNPLAN.json`. The portable Forward entry is:

```bash
python paper_cases/dna/DNA01/models/nt-v2/forward/scripts/encode_dna01_transformer_models.py \
  --model nt --strict-repeat \
  --object-table constraints/dna/DNA01/refset/DNA01_OBJECT_TABLE.tsv \
  --fasta constraints/dna/DNA01/refset/records.fasta \
  --output work/DNA01/nt-v2/forward
```

Replace angle-bracketed model locations with the exact resources identified under `paper_cases/model_resources/dna/nt-v2/`. Do not change the RefSet population or extraction semantics.
