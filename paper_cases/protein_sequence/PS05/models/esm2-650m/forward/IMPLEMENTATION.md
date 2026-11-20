# PS05 x ESM-2 650M Forward

This case regenerates the PS05 scorer input from the released 321-object FASTA and the registered ESM-2 checkpoint. It does not implement PS05 scoring.

- Public Forward adapter: [`../../../shared_forward/ps05_forward_public.py`](../../../shared_forward/ps05_forward_public.py)
- Public scoring command: [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)
- Constraint contract and scorer: [`constraints/protein_sequence/PS05/README.md`](../../../../../../constraints/protein_sequence/PS05/README.md)

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS05/shared_forward/ps05_forward_public.py \
  --model esm2 \
  --fasta constraints/protein_sequence/PS05/refset/PS05_ENDPOINT_SEQUENCES.faa \
  --checkpoint /path/to/esm2_t33_650M_UR50D.pt \
  --output-dir work/PS05/esm2-650m/forward
```

The adapter exports `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`. It uses ESM-2 layer 33 and residue-mean pooling, matching the registered paper interface. Continue with the public scoring command.
