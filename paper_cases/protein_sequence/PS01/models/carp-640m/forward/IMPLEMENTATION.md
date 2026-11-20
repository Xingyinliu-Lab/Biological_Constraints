# carp-640m Forward implementation

The current portable entry point is [`../../../shared_forward/ps01_forward_public.py`](../../../shared_forward/ps01_forward_public.py). The model-local script below is a thin wrapper around the same parameterized adapter; it contains no server-specific input, checkpoint or work path.

- `ps01_carp640_forward.py`+

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS01/shared_forward/ps01_forward_public.py \
  --model carp \
  --fasta constraints/protein_sequence/PS01/refset/FORMAL_SEQUENCES.faa.gz \
  --checkpoint /path/to/TORCH_HOME/hub/checkpoints/carp_640M.pt \
  --output-dir work/PS01/carp-640m/forward
```

The command exports `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json` for all 2,316 released sequences.

The validated pooling expression is `h[j, :len(seqs[i])].float().mean(0)`. The 2026-10-19 correction replaced an erroneous occurrence-ID-length slice; the affected CARP-640M case was rerun and rescored.
