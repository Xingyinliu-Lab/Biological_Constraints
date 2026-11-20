
# DNA01 × caduceus-ps: execution commands

The registered execution parameters are in `RUNPLAN.json`. The portable Forward entry is:

```bash
python paper_cases/dna/DNA01/models/caduceus-ps/forward/scripts/encode_dna01_transformer_models.py \
  --model caduceus_ps --strict-repeat \
  --object-table constraints/dna/DNA01/refset/DNA01_OBJECT_TABLE.tsv \
  --fasta constraints/dna/DNA01/refset/records.fasta \
  --output work/DNA01/caduceus-ps/forward
```

Replace angle-bracketed model locations with the exact resources identified under `paper_cases/model_resources/dna/caduceus-ps/`. Do not change the RefSet population or extraction semantics.
