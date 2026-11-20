# ESM-C 600M Forward for PS02A

The executed implementation is [`scripts/ps02a_forward.py`](scripts/ps02a_forward.py).
It reads all 13,720 released Forward FASTA records, excludes special/padding tokens,
mean-pools the final sequence hidden states and writes
`embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`.

For a current run, use:

```bash
# Run from the repository root.
python paper_cases/protein_sequence/PS02A/models/esmc-600m/forward/scripts/ps02a_forward.py esmc --fasta constraints/protein_sequence/PS02A/refset/endpoint_sequences.faa --checkpoint /models/esmc-600m --workdir work/PS02A/esmc-600m/forward
```

This writes `work/PS02A/esmc-600m/forward/esmc/embeddings.safetensors` and its index. Retain the frozen source revision, pooling, full coverage and output schema.

Pass the saved tensors and index to
[`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
