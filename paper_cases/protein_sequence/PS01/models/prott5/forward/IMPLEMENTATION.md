# prott5 Forward implementation

The current portable entry point is [`../../../shared_forward/ps01_forward_public.py`](../../../shared_forward/ps01_forward_public.py). The model-local script below is a thin wrapper around the same parameterized adapter; it contains no server-specific input, checkpoint or work path.

- `ps01_forward.py`
- `run_forward.sh`+

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS01/shared_forward/ps01_forward_public.py \
  --model prott5 \
  --fasta constraints/protein_sequence/PS01/refset/FORMAL_SEQUENCES.faa.gz \
  --checkpoint /path/to/prot_t5_xl_half_uniref50-enc \
  --output-dir work/PS01/prott5/forward
```

The command exports `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json` for all 2,316 released sequences.
