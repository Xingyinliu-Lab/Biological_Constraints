# CARP-640M Forward for PS02A

The executed implementation is
[`scripts/ps02a_carp640_forward.py`](scripts/ps02a_carp640_forward.py). It reads
all 13,720 released Forward sequences, obtains the final ByteNet representation, excludes
padding and mean-pools residues. Outputs are `embeddings.safetensors`,
`embedding_index.tsv` and `FORWARD_REPORT.json`.

For a current run, provide the released FASTA, a standard
`TORCH_HOME/hub/checkpoints/carp_640M.pt` checkpoint and a new work directory:

```bash
# Run from the repository root.
python paper_cases/protein_sequence/PS02A/models/carp-640m/forward/scripts/ps02a_carp640_forward.py --fasta constraints/protein_sequence/PS02A/refset/endpoint_sequences.faa --checkpoint /torch-cache/hub/checkpoints/carp_640M.pt --workdir work/PS02A/carp-640m/forward
```

The script derives `TORCH_HOME` from the checkpoint path. Model identity and environment are documented in
[`paper_cases/model_resources/protein_sequence/carp-640m/`](../../../../../model_resources/protein_sequence/carp-640m/README.md).
Do not alter the model layer, pooling, object population or schema.

The current implementation pools each sequence through its actual amino-acid
length (`len(seqs[i])`), excluding padding.

Pass `work/PS02A/carp-640m/forward/carp/embeddings.safetensors` and its index to
[`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
