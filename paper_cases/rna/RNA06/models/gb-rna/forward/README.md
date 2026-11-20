# RNA06 Forward: gb-rna

Run from the repository root. The encoder and exporter use one shared work directory:

```bash
mkdir -p work/RNA06/gb-rna/forward
python paper_cases/rna/RNA06/models/gb-rna/forward/encode_gb_nucleotide.py --fasta constraints/rna/RNA06/refset/rna_sequences.fasta --model-dir <MODEL_DIR> --output work/RNA06/gb-rna/forward/SEQUENCE_MATRICES.pt
python paper_cases/rna/RNA06/shared_forward/export_formal_positions.py --matrices work/RNA06/gb-rna/forward/SEQUENCE_MATRICES.pt --output work/RNA06/gb-rna/forward/REPRESENTATIONS.npz
```

The first command performs model-specific Forward; the shared exporter only selects the frozen biological positions and writes `sequence_sha256 + position0 + vectors`.

The large representation file is not redistributed. Its size and SHA256 are in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json). Model environment and checkpoint identity are maintained under `paper_cases/model_resources/rna/gb-rna/`.
