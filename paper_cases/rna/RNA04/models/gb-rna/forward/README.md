# RNA04 × gb-rna Forward

This Forward uses the released endpoint FASTA directly and writes the `sequence_sha256` plus `vectors` NPZ required by the [public Constraint](../../../../../../constraints/rna/RNA04/README.md). Run from the repository root after obtaining the checkpoint described in [`paper_cases/model_resources/rna/gb-rna/`](../../../../../model_resources/rna/gb-rna/):

```bash
python paper_cases/rna/RNA04/models/gb-rna/forward/encode.py \
  --fasta constraints/rna/RNA04/refset/RNA04_unique_endpoints.fasta \
  --model-dir /path/to/GB-RNA-650M \
  --output work/RNA04/gb-rna/EMBEDDINGS.npz
```

The adapter pools valid nucleotide states only, excludes special tokens, and averages complete 1,022-nt windows to a whole-CDS representation. Score the resulting NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).

The current Forward implementation requires every 1,022-nt-or-shorter window to produce exactly one valid nucleotide token per input nucleotide before pooling.
