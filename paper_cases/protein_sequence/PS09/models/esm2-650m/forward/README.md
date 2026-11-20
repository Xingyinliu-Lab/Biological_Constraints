# PS09 x ESM-2 650M Forward

The executed Forward generates one residue representation for every formal PS09 endpoint. The output is `endpoint_embeddings.safetensors`, keyed by the 15,708 endpoint IDs defined by the Constraint.

Current configurable Forward source: [`../../../shared_forward/ps09_four_model_forward.py`](../../../shared_forward/ps09_four_model_forward.py).

```bash
python paper_cases/protein_sequence/PS09/shared_forward/ps09_four_model_forward.py \
  --model esm2-650m \
  --checkpoint /path/to/esm2_t33_650M_UR50D.pt \
  --refset constraints/protein_sequence/PS09/refset \
  --output-dir work/PS09/esm2-650m/forward
```

Use the registered model environment and checkpoint documented by this case. If fair-esm is not installed in that environment, pass `--esm2-import-path` to its Python package directory. The entry does not require the author's server layout. After Forward, use the public command in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
