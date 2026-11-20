# PS05 x ESM3 Sequence Forward

This case preserves the executed PS05 model-forward implementation. It generates the representation tensor and ID index accepted by the PS05 public scorer; it does not implement PS05 scoring.

- Forward code: [`scripts/esm3_seq_mean_extract.py`](scripts/esm3_seq_mean_extract.py)
- Public scoring command: [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)
- Constraint contract and scorer: [`constraints/protein_sequence/PS05/README.md`](../../../../../../constraints/protein_sequence/PS05/README.md)

Run from the repository root with the registered ESM3 sequence-track checkpoint available in the documented model environment:

```bash
python paper_cases/protein_sequence/PS05/models/esm3-seq/forward/scripts/esm3_seq_mean_extract.py \
  --fasta constraints/protein_sequence/PS05/refset/PS05_ENDPOINT_SEQUENCES.faa \
  --output-dir work/PS05/esm3-seq/forward \
  --id-column domain_id
```

The Forward output must preserve the released PS05 endpoint IDs exactly. Checkpoint identity is described by the shared model-resource records.
