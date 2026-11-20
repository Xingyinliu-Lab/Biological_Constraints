# esm3-seq Forward implementation

The retained scripts are the actual paper Forward implementation for the 2,316 released PS01 sequences. They produce whole-protein vectors and the safetensors/index pair accepted by the public scorer.

- `esm3_seq_mean_extract.py`

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS01/models/esm3-seq/forward/scripts/esm3_seq_mean_extract.py \
  --fasta constraints/protein_sequence/PS01/refset/FORMAL_SEQUENCES.faa.gz \
  --output-dir work/PS01/esm3-seq/forward \
  --id-column occurrence_id
```

This produces `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json` for all 2,316 released sequences.
