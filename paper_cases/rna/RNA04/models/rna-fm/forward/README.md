# RNA04 × rna-fm Forward

This Forward uses the released endpoint FASTA directly and writes the `sequence_sha256` plus `vectors` NPZ required by the [public Constraint](../../../../../../constraints/rna/RNA04/README.md). Run from the repository root after obtaining the checkpoint and source checkout described in [`paper_cases/model_resources/rna/rna-fm/`](../../../../../model_resources/rna/rna-fm/):

```bash
python paper_cases/rna/RNA04/models/rna-fm/forward/encode.py \
  --fasta constraints/rna/RNA04/refset/RNA04_unique_endpoints.fasta \
  --weights /path/to/RNA-FM_pretrained.pth \
  --source-root /path/to/RNA-FM \
  --output work/RNA04/rna-fm/EMBEDDINGS.npz
```

The adapter uses overlapping 1,022-nt windows (stride 766), averages hidden states by nucleotide coordinate, excludes BOS/EOS, then takes an equal-nucleotide whole-CDS mean. Score the resulting NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
