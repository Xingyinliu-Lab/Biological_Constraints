# ESM-2 Forward for PS02A

The executed implementation is [`scripts/ps02a_forward.py`](scripts/ps02a_forward.py).
It reads all 13,720 released Forward FASTA records, excludes special tokens, mean-pools
residue representations from layer 33 and writes:

- `embeddings.safetensors`: one float32 vector tensor per object;
- `embedding_index.tsv`: `occurrence_id`, `embedding_key`, sequence length,
  embedding dimension and vector SHA256;
- `FORWARD_REPORT.json`: coverage, checkpoint identity and output identity.

For a current run, use the released FASTA directly:

```bash
# Run from the repository root.
python paper_cases/protein_sequence/PS02A/models/esm2-650m/forward/scripts/ps02a_forward.py esm2 --fasta constraints/protein_sequence/PS02A/refset/endpoint_sequences.faa --checkpoint /models/esm2_t33_650M_UR50D.pt --workdir work/PS02A/esm2-650m/forward
```

This writes `work/PS02A/esm2-650m/forward/esm2/embeddings.safetensors` and its index. Do not alter layer, pooling, object membership or output schema.

Pass the safetensors file and index to the command in
[`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
