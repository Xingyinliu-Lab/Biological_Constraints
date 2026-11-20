# Forward implementation

Run from the repository root. Use the current parameterized entry point:

```bash
python paper_cases/protein_sequence/PS08/shared_forward/ps08_forward_public.py --model esm2 --fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz --checkpoint /models/esm2_t33_650M_UR50D.pt --output-dir work/PS08/esm2-650m
```

It mean-pools ESM-2 layer 33 over residues, excluding special tokens, and writes the safetensors/index pair. Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
