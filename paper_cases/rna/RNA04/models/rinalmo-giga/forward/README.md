# RNA04 × rinalmo-giga Forward

This Forward uses the released endpoint FASTA directly and writes the `sequence_sha256` plus `vectors` NPZ required by the [public Constraint](../../../../../../constraints/rna/RNA04/README.md). Run from the repository root after obtaining the source checkout described in [`paper_cases/model_resources/rna/rinalmo-giga/`](../../../../../model_resources/rna/rinalmo-giga/):

```bash
python paper_cases/rna/RNA04/models/rinalmo-giga/forward/encode.py \
  --fasta constraints/rna/RNA04/refset/RNA04_unique_endpoints.fasta \
  --repo /path/to/rinalmo-source \
  --out work/RNA04/rinalmo-giga/EMBEDDINGS.npz
```

The adapter uses overlapping 1,022-nt windows (stride 766), averages hidden states by nucleotide coordinate, excludes CLS/EOS, then takes an equal-nucleotide whole-CDS mean. Score the resulting NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
