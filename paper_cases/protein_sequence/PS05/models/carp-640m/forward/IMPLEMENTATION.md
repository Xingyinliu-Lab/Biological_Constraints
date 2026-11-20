# PS05 x CARP-640M Forward

This case regenerates the PS05 scorer input from the released 321-object FASTA and the registered CARP-640M checkpoint. It does not implement PS05 scoring.

- Public Forward adapter: [`../../../shared_forward/ps05_forward_public.py`](../../../shared_forward/ps05_forward_public.py)
- Public scoring command: [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)
- Constraint contract and scorer: [`constraints/protein_sequence/PS05/README.md`](../../../../../../constraints/protein_sequence/PS05/README.md)

Run from the repository root. The checkpoint must retain the standard `<TORCH_HOME>/hub/checkpoints/carp_640M.pt` layout used by the CARP loader.

```bash
python paper_cases/protein_sequence/PS05/shared_forward/ps05_forward_public.py \
  --model carp \
  --fasta constraints/protein_sequence/PS05/refset/PS05_ENDPOINT_SEQUENCES.faa \
  --checkpoint /path/to/TORCH_HOME/hub/checkpoints/carp_640M.pt \
  --output-dir work/PS05/carp-640m/forward
```

The adapter extracts the final ByteNet residue states, applies residue-mean pooling and exports `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json`. Continue with the public scoring command.
