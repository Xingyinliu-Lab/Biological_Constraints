# PS09 x ESM3 sequence Forward

The executed Forward generates one residue representation for every formal PS09 endpoint. The output is `endpoint_embeddings.safetensors`, keyed by the 15,708 endpoint IDs defined by the Constraint. It does not use global sequence mean pooling: each submitted vector is the hidden state for its specific frozen residue coordinate.

Executed Forward source: [`scripts/esm3_seq_residue_extract.py`](scripts/esm3_seq_residue_extract.py).

```bash
python paper_cases/protein_sequence/PS09/models/esm3-seq/forward/scripts/esm3_seq_residue_extract.py --membership constraints/protein_sequence/PS09/refset/FORMAL_BLOCK_MEMBERSHIP.tsv.gz --fasta constraints/protein_sequence/PS09/refset/SELECTED_PROTEINS.faa --output-dir work/PS09/esm3-seq/forward
```

Do not change token-to-residue mapping, windowing, layer selection or residue-level readout. After Forward, use the public command in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
