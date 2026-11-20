# RNA10 Forward: orthrus-4-track

The recovered Python implementation encodes every released WT and mutant sequence and writes the sequence-level NPZ required by the public RNA10 measurement builder. The generated NPZ is identified in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json) and is not redistributed.

Run from the repository root using the registered [Orthrus model resource](../../../../../model_resources/rna/orthrus-4-track/README.md), [environment](../../../../../model_resources/rna/orthrus-4-track/ENVIRONMENT.md), and [checkpoint identity](../../../../../model_resources/rna/orthrus-4-track/checkpoint_identity.json):

```bash
python paper_cases/rna/RNA10/models/orthrus-4-track/forward/encode_orthrus_4track.py \
  --fasta constraints/rna/RNA10/refset/RNA10_UNIQUE_INPUTS.fasta \
  --output work/RNA10/orthrus-4-track/forward/REPRESENTATIONS.npz \
  --model-id antichronology/orthrus-4-track \
  --device cuda
```

This RNA10 implementation exposes only the formal whole-sequence representation. Continue with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
