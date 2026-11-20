# ESM3 sequence-track Forward for PS02A

Run the executed extractor in the ESM3 model environment:

```bash
python paper_cases/protein_sequence/PS02A/models/esm3-seq/forward/scripts/esm3_seq_mean_extract.py \
  --fasta constraints/protein_sequence/PS02A/refset/endpoint_sequences.faa \
  --output-dir work/PS02A/esm3-seq/forward \
  --id-column occurrence_id
```

It uses the frozen sequence-only track, bfloat16 compute, 1,023-residue windows,
767-residue stride (256-residue overlap), tail-window alignment and residue-wise
mean pooling across overlapping windows. It writes
`embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json` for
all 13,720 released Forward objects. The index maps formal object IDs to tensor keys.

Pass those two representation files to
[`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
