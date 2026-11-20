# RNA04 × orthrus-4-track Forward

This Forward uses the released endpoint FASTA directly and writes the `sequence_sha256` plus `vectors` NPZ required by the [public Constraint](../../../../../../constraints/rna/RNA04/README.md). Run from the repository root after obtaining the checkpoint described in [`paper_cases/model_resources/rna/orthrus-4-track/`](../../../../../model_resources/rna/orthrus-4-track/):

```bash
python paper_cases/rna/RNA04/models/orthrus-4-track/forward/encode.py \
  --fasta constraints/rna/RNA04/refset/RNA04_unique_endpoints.fasta \
  --mode global \
  --model-id /path/or/huggingface-id/of/orthrus-4-track \
  --output work/RNA04/orthrus-4-track/EMBEDDINGS.npz
```

`global` requests the model's whole-sequence representation after the adapter validates the SHA-addressed A/C/G/U input. Score the resulting NPZ with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
