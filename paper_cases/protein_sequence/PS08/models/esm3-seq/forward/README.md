# Forward implementation

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS08/models/esm3-seq/forward/esm3_seq_mean_extract.py --fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz --id-column endpoint_id --output-dir work/PS08/esm3-seq
```

This sequence-only extractor excludes special tokens, applies its documented window/overlap reconstruction, and residue-mean-pools the result. Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
