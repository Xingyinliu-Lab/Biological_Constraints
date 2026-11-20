
# DNA01 × dnabert-2: execution commands

The registered execution parameters are in `RUNPLAN.json`. The portable Forward entry is the registered Hugging Face Direct implementation (`zhihan1996/DNABERT-2-117M`, revision `7bce263b15377fc15361f52cfab88f8b586abda0`, float32 and standard-attention fallback):

```bash
python paper_cases/dna/DNA01/models/dnabert-2/forward/scripts/encode_dna01_transformer_models.py \
  --model dnabert2 --strict-repeat \
  --dnabert2-model-path <DNABERT2_SNAPSHOT_DIR> \
  --object-table constraints/dna/DNA01/refset/DNA01_OBJECT_TABLE.tsv \
  --fasta constraints/dna/DNA01/refset/records.fasta \
  --output work/DNA01/dnabert-2/forward
```

Replace angle-bracketed model locations with the exact resources identified under `paper_cases/model_resources/dna/dnabert-2/`. Do not change the RefSet population or extraction semantics.
