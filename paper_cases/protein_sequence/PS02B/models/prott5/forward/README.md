# ProtT5 Forward for PS02B

Run from the repository root using the registered model checkpoint:

```bash
python paper_cases/protein_sequence/PS02B/models/prott5/forward/scripts/ps02b_forward.py prott5 --fasta constraints/protein_sequence/PS02B/refset/endpoint_sequences.faa --checkpoint /path/to/ProtT5 --workdir work/PS02B/prott5/forward
```

The executed implementation is [`scripts/ps02b_forward.py`](scripts/ps02b_forward.py). It reads all 3,937 released FASTA objects, applies the registered ProtT5 tokenization/extraction and mean pooling, and writes `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`.

The command above sets the released input and output paths. Use the public scorer in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
