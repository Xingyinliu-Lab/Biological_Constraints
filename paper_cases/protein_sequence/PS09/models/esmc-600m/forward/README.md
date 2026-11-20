# PS09 x ESM-C 600M Forward

The executed Forward generates one residue representation for every formal PS09 endpoint. The output is `endpoint_embeddings.safetensors`, keyed by the 15,708 endpoint IDs defined by the Constraint.

Current configurable Forward source: [`../../../shared_forward/ps09_four_model_forward.py`](../../../shared_forward/ps09_four_model_forward.py).

```bash
python paper_cases/protein_sequence/PS09/shared_forward/ps09_four_model_forward.py \
  --model esmc-600m \
  --checkpoint /path/to/biohub_ESMC-600M \
  --refset constraints/protein_sequence/PS09/refset \
  --output-dir work/PS09/esmc-600m/forward
```

Use the registered model environment and checkpoint documented by this case. The entry accepts explicit paths and does not require the author's server layout. After Forward, use the public command in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
