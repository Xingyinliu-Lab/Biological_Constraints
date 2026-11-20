# PS05 x ESM-C 600M Forward

This case regenerates the PS05 scorer input from the released 321-object FASTA and the registered ESM-C checkpoint directory. It does not implement PS05 scoring.

- Public Forward adapter: [`../../../shared_forward/ps05_forward_public.py`](../../../shared_forward/ps05_forward_public.py)
- Public scoring command: [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)
- Constraint contract and scorer: [`constraints/protein_sequence/PS05/README.md`](../../../../../../constraints/protein_sequence/PS05/README.md)

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS05/shared_forward/ps05_forward_public.py \
  --model esmc \
  --fasta constraints/protein_sequence/PS05/refset/PS05_ENDPOINT_SEQUENCES.faa \
  --checkpoint /path/to/ESMC-600M \
  --output-dir work/PS05/esmc-600m/forward
```

The adapter excludes special-token positions, applies residue-mean pooling and exports `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`. Continue with the public scoring command.
