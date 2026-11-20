# RNA10 Forward: gb-rna

The recovered Python implementation encodes every released WT and mutant sequence and writes the sequence-level NPZ required by the public RNA10 measurement builder. The generated NPZ is identified in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json) and is not redistributed.

Run from the repository root using the registered [GB-RNA model resource](../../../../../model_resources/rna/gb-rna/README.md), [environment](../../../../../model_resources/rna/gb-rna/ENVIRONMENT.md), and [checkpoint identity](../../../../../model_resources/rna/gb-rna/checkpoint_identity.json):

```bash
python paper_cases/rna/RNA10/models/gb-rna/forward/encode_gb_rna_global.py \
  --fasta constraints/rna/RNA10/refset/RNA10_UNIQUE_INPUTS.fasta \
  --model-dir <GB_RNA_MODEL_DIRECTORY> \
  --output work/RNA10/gb-rna/forward/REPRESENTATIONS.npz \
  --device cuda
```

Do not use `--limit` for a formal run. Continue with the commands in [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
