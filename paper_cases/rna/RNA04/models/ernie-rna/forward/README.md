# RNA04 × ernie-rna Forward

This Forward uses the released endpoint FASTA directly and writes the `sequence_sha256` plus `vectors` NPZ required by the [public Constraint](../../../../../../constraints/rna/RNA04/README.md). Run from the repository root after obtaining the checkpoint described in [`paper_cases/model_resources/rna/ernie-rna/`](../../../../../model_resources/rna/ernie-rna/):

```bash
python paper_cases/rna/RNA04/models/ernie-rna/forward/encode.py \
  --fasta constraints/rna/RNA04/refset/RNA04_unique_endpoints.fasta \
  --model-dir /path/to/ERNIE-RNA-HF \
  --output work/RNA04/ernie-rna/EMBEDDINGS.npz
```

The adapter uses overlapping 1,022-nt windows (stride 766), averages hidden states by nucleotide coordinate across overlaps, excludes special tokens, then takes an equal-nucleotide whole-CDS mean. It verifies one nucleotide token per input nucleotide. Score the resulting NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
