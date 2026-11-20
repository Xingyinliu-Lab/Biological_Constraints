# RNA10 Forward: rna-fm

The recovered Python implementation encodes every released WT and mutant sequence and writes the sequence-level NPZ required by the public RNA10 measurement builder. The generated NPZ is identified in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json) and is not redistributed.

Run from the repository root using the registered [RNA-FM model resource](../../../../../model_resources/rna/rna-fm/README.md), [environment](../../../../../model_resources/rna/rna-fm/ENVIRONMENT.md), and [checkpoint identity](../../../../../model_resources/rna/rna-fm/checkpoint_identity.json):

```bash
python paper_cases/rna/RNA10/models/rna-fm/forward/encode_rnafm_global.py \
  --fasta constraints/rna/RNA10/refset/RNA10_UNIQUE_INPUTS.fasta \
  --weights <RNA_FM_WEIGHTS> \
  --source-root <RNA_FM_SOURCE_ROOT> \
  --output work/RNA10/rna-fm/forward/REPRESENTATIONS.npz \
  --device cuda
```

Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
