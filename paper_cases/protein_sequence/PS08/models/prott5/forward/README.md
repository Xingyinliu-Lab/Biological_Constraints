# Forward implementation

Run from the repository root. Use the current parameterized entry point:

```bash
python paper_cases/protein_sequence/PS08/shared_forward/ps08_forward_public.py --model prott5 --fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz --checkpoint /models/prot_t5_xl_half_uniref50-enc --output-dir work/PS08/prott5
```

It records a deterministic combined SHA256 for the supplied checkpoint directory, maps `U/O` to `X`, excludes special tokens, and residue-mean-pools the encoder output. Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
