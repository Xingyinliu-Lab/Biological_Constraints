# ESM-2 Forward for PS02B

Run from the repository root using the registered model checkpoint:

```bash
python paper_cases/protein_sequence/PS02B/models/esm2-650m/forward/scripts/ps02b_forward.py esm2 --fasta constraints/protein_sequence/PS02B/refset/endpoint_sequences.faa --checkpoint /path/to/esm2_t33_650M_UR50D.pt --workdir work/PS02B/esm2-650m/forward
```

The executed implementation is [`scripts/ps02b_forward.py`](scripts/ps02b_forward.py). It reads all 3,937 released FASTA objects, extracts layer-33 residue representations, excludes special tokens, mean-pools one vector per domain and writes `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`.

The command above sets the released input and output paths. Use the public scorer in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
