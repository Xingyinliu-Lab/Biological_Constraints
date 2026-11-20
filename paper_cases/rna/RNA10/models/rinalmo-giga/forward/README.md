# RNA10 Forward: rinalmo-giga

The recovered Python implementation encodes every released WT and mutant sequence and writes the sequence-level NPZ required by the public RNA10 measurement builder. The generated NPZ is identified in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json) and is not redistributed.

Run from the repository root using the registered [RiNALMo model resource](../../../../../model_resources/rna/rinalmo-giga/README.md), [environment](../../../../../model_resources/rna/rinalmo-giga/ENVIRONMENT.md), and [checkpoint identity](../../../../../model_resources/rna/rinalmo-giga/checkpoint_identity.json):

```bash
python paper_cases/rna/RNA10/models/rinalmo-giga/forward/encode_official_rinalmo_global.py \
  --fasta constraints/rna/RNA10/refset/RNA10_UNIQUE_INPUTS.fasta \
  --repo <RINALMO_REPOSITORY> \
  --out work/RNA10/rinalmo-giga/forward/REPRESENTATIONS.npz \
  --device cuda
```

Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
