# RNA06 Forward: rna-fm

Run from the repository root. The encoder and exporter use one shared work directory:

```bash
mkdir -p work/RNA06/rna-fm/forward
python paper_cases/rna/RNA06/models/rna-fm/forward/encode_rnafm_nucleotide.py --fasta constraints/rna/RNA06/refset/rna_sequences.fasta --weights <WEIGHTS> --source-root <RNA_FM_SOURCE> --output work/RNA06/rna-fm/forward/SEQUENCE_MATRICES.pt
python paper_cases/rna/RNA06/shared_forward/export_formal_positions.py --matrices work/RNA06/rna-fm/forward/SEQUENCE_MATRICES.pt --output work/RNA06/rna-fm/forward/REPRESENTATIONS.npz
```

The first command performs model-specific Forward; the shared exporter only selects the frozen biological positions and writes `sequence_sha256 + position0 + vectors`.

The large representation file is not redistributed. Its size and SHA256 are in [`../scoring_call/SCORER_INPUT_BINDING.json`](../scoring_call/SCORER_INPUT_BINDING.json). Model environment and checkpoint identity are maintained under `paper_cases/model_resources/rna/rna-fm/`.
