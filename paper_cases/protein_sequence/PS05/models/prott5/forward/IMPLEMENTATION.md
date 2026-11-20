# PS05 x ProtT5 Forward

This case regenerates the PS05 scorer input from the released 321-object FASTA and the registered ProtT5 checkpoint. It does not implement PS05 scoring.

- Public Forward adapter: [`../../../shared_forward/ps05_forward_public.py`](../../../shared_forward/ps05_forward_public.py)
- Public scoring command: [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)
- Constraint contract and scorer: [`constraints/protein_sequence/PS05/README.md`](../../../../../../constraints/protein_sequence/PS05/README.md)

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS05/shared_forward/ps05_forward_public.py \
  --model prott5 \
  --fasta constraints/protein_sequence/PS05/refset/PS05_ENDPOINT_SEQUENCES.faa \
  --checkpoint /path/to/prot_t5_xl_half_uniref50-enc \
  --output-dir work/PS05/prott5/forward
```

The adapter applies the registered residue tokenization, maps `U`, `Z`, `O` and `B` to `X`, removes special-token positions from pooling and exports `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`. Continue with the public scoring command.
