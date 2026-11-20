# RNA05 Forward: orthrus-4-track

`encode.py` is the current verified paper Forward implementation. From the repository root:

```bash
mkdir -p work/RNA05/orthrus-4-track/forward
python paper_cases/rna/RNA05/models/orthrus-4-track/forward/encode.py \
  --fasta constraints/rna/RNA05/refset/RNA05_06_unique_inputs.fasta \
  --mode positions \
  --output work/RNA05/orthrus-4-track/forward/REPRESENTATIONS.npz
```

`--mode positions` is required for the RNA05 nucleotide-level contract. The output is keyed by sequence SHA256 with one `(sequence_length, d)` array per sequence.

The large representation file is not redistributed here. Its byte size and SHA256 are recorded in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json). Model environment and checkpoint identity are maintained under `paper_cases/model_resources/rna/orthrus-4-track/`.

After Forward, use [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
