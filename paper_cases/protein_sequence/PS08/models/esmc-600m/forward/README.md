# Forward implementation

Run from the repository root. Use the current parameterized entry point:

```bash
python paper_cases/protein_sequence/PS08/shared_forward/ps08_forward_public.py --model esmc --fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz --checkpoint /models/biohub_ESMC-600M --output-dir work/PS08/esmc-600m
```

The script records a deterministic combined SHA256 for the supplied checkpoint directory, excludes special tokens before residue-mean pooling, and then continues with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
