# esmc-600m Forward implementation

The current portable entry point is [`../../../shared_forward/ps01_forward_public.py`](../../../shared_forward/ps01_forward_public.py). The model-local script below is a thin wrapper around the same parameterized adapter; it contains no server-specific input, checkpoint or work path.

- `ps01_forward.py`
- `run_forward.sh`+

Run from the repository root:

```bash
python paper_cases/protein_sequence/PS01/shared_forward/ps01_forward_public.py \
  --model esmc \
  --fasta constraints/protein_sequence/PS01/refset/FORMAL_SEQUENCES.faa.gz \
  --checkpoint /path/to/ESMC-600M \
  --output-dir work/PS01/esmc-600m/forward
```

The command exports `embeddings.safetensors`, `embedding_index.tsv` and `FORWARD_REPORT.json` for all 2,316 released sequences.
