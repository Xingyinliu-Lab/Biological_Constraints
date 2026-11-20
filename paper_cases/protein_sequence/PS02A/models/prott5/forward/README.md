# ProtT5 Forward for PS02A

The executed implementation is [`scripts/ps02a_forward.py`](scripts/ps02a_forward.py).
It reads all 13,720 released Forward FASTA records, applies the frozen `U/O → X`
normalization, excludes padding, mean-pools final hidden states and writes
`embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`.

For a current run, use:

```bash
# Run from the repository root.
python paper_cases/protein_sequence/PS02A/models/prott5/forward/scripts/ps02a_forward.py prott5 --fasta constraints/protein_sequence/PS02A/refset/endpoint_sequences.faa --checkpoint /models/prott5_direct --workdir work/PS02A/prott5/forward
```

This writes `work/PS02A/prott5/forward/prott5/embeddings.safetensors` and its index. Keep population, normalization, pooling and output schema unchanged.

The index maps every `occurrence_id` to its tensor key. Pass both files to
[`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
