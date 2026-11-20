# ESM-C 600M Forward for PS02B

Run from the repository root using the registered model checkpoint:

```bash
python paper_cases/protein_sequence/PS02B/models/esmc-600m/forward/scripts/ps02b_forward.py esmc --fasta constraints/protein_sequence/PS02B/refset/endpoint_sequences.faa --checkpoint /path/to/ESMC-600M --workdir work/PS02B/esmc-600m/forward
```

The executed implementation is [`scripts/ps02b_forward.py`](scripts/ps02b_forward.py). It reads all 3,937 released FASTA objects, performs the registered ESM-C extraction and mean pooling, and writes `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`.

The command above sets the released input and output paths. Use the public scorer in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
