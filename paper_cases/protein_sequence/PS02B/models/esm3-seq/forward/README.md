# ESM3 sequence-track Forward for PS02B

Run from the repository root using the registered model checkpoint:

```bash
python paper_cases/protein_sequence/PS02B/models/esm3-seq/forward/scripts/esm3_seq_mean_extract.py --fasta constraints/protein_sequence/PS02B/refset/endpoint_sequences.faa --id-column domain_id --output-dir work/PS02B/esm3-seq/forward
```

The executed implementation is [`scripts/esm3_seq_mean_extract.py`](scripts/esm3_seq_mean_extract.py). It uses the registered ESM3 sequence-track checkpoint and extraction settings to emit one whole-domain vector for every released FASTA object, with an index and Forward identity report.

The command above sets the released input and output paths. Use the public scorer in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
