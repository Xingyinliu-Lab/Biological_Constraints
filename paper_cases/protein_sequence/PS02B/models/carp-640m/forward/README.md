# CARP-640M Forward for PS02B

Run from the repository root using the registered model checkpoint:

```bash
python paper_cases/protein_sequence/PS02B/models/carp-640m/forward/scripts/ps02b_carp640_forward.py --fasta constraints/protein_sequence/PS02B/refset/endpoint_sequences.faa --checkpoint /path/to/TORCH_HOME/hub/checkpoints/carp_640M.pt --workdir work/PS02B/carp-640m/forward
```

The executed implementation is [`scripts/ps02b_carp640_forward.py`](scripts/ps02b_carp640_forward.py). It reads all 3,937 released FASTA objects, applies the registered CARP extraction, mean-pools residues over each object's actual sequence length, and writes one indexed vector per domain plus a Forward report.

The command above sets the released input and output paths. Use the public scorer in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
