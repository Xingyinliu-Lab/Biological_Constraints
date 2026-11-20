# RNA06 Forward: orthrus-4-track

Run from the repository root. The encoder and exporter use one shared work directory:

```bash
mkdir -p work/RNA06/orthrus-4-track/forward
python paper_cases/rna/RNA06/models/orthrus-4-track/forward/encode_orthrus_4track.py --fasta constraints/rna/RNA06/refset/rna_sequences.fasta --mode positions --output work/RNA06/orthrus-4-track/forward/SEQUENCE_MATRICES.npz
python paper_cases/rna/RNA06/shared_forward/export_formal_positions.py --matrices work/RNA06/orthrus-4-track/forward/SEQUENCE_MATRICES.npz --output work/RNA06/orthrus-4-track/forward/REPRESENTATIONS.npz
```

The first command performs model-specific Forward; the shared exporter only selects the frozen biological positions and writes `sequence_sha256 + position0 + vectors`.

The large representation file is not redistributed. Its size and SHA256 are in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json). Model environment and checkpoint identity are maintained under `paper_cases/model_resources/rna/orthrus-4-track/`.
