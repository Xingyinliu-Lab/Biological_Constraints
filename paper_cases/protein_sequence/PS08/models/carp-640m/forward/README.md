# Forward implementation

Run from the repository root. Use the current parameterized entry point:

```bash
python paper_cases/protein_sequence/PS08/shared_forward/ps08_forward_public.py --model carp --fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz --checkpoint /models/carp_640M.pt --output-dir work/PS08/carp-640m
```

The script loads the specified checkpoint file directly and residue-mean-pools the final layer. Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
