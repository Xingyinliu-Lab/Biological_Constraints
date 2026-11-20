# RNA05 Forward: ernie-rna

`encode.py` is the current verified paper Forward implementation. From the repository root:

```bash
mkdir -p work/RNA05/ernie-rna/forward
python paper_cases/rna/RNA05/models/ernie-rna/forward/encode.py \
  --fasta constraints/rna/RNA05/refset/RNA05_06_unique_inputs.fasta \
  --model-dir <ERNIE_RNA_CHECKPOINT> \
  --output work/RNA05/ernie-rna/forward/REPRESENTATIONS.npz
```

It writes a NumPy `REPRESENTATIONS.npz`, keyed by sequence SHA256 with one `(sequence_length, d)` array per sequence.

The large representation file is not redistributed here. Its byte size and SHA256 are recorded in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json). Model environment and checkpoint identity are maintained under `paper_cases/model_resources/rna/ernie-rna/`.

After Forward, use [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
