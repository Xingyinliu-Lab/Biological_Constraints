# RNA01-02 × RNA-FM Forward

This directory preserves the Forward implementation used to produce the registered RNA01-02 representations. It does not contain a scoring algorithm.

## Inputs and outputs

The Forward input and representation contract are defined by the [RNA01-02 Constraint](../../../../../../constraints/rna/RNA01-02/README.md). The implementation must export the full L16K vectors and two calibration-repeat files with the exact object IDs required by that contract.

## Implementation

- Executed FULL-package entry: [`executed_full/encode.py`](executed_full/encode.py)
- Recovered supporting scripts: [`scripts/`](scripts/)
- Model identity and environment: [`paper_cases/model_resources/rna/rna-fm/`](../../../../../model_resources/rna/rna-fm/)
- Public scorer invocation: [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)

The public Constraint scorer performs calibration, relation measurement, aggregation and confidence-interval calculation. Result-local scoring implementations are not part of the release.

## Current runnable sequence

From the repository root, first run the RefSet materialization command in the Constraint README. Then run the released RNA-FM extractor once for L16K and twice for the calibration FASTA:

```bash
mkdir -p work/RNA01-02/rna-fm
python paper_cases/rna/RNA01-02/models/rna-fm/forward/executed_full/encode.py --fasta work/RNA01-02/forward_input/L16K_UNIQUE_SEQUENCES.fasta --weights /path/to/rnafm_weights.pt --source-root /path/to/RNA-FM --output work/RNA01-02/rna-fm/L16K_VECTORS.npz
python paper_cases/rna/RNA01-02/models/rna-fm/forward/executed_full/encode.py --fasta work/RNA01-02/forward_input/CALIBRATION_UNIQUE_SEQUENCES.fasta --weights /path/to/rnafm_weights.pt --source-root /path/to/RNA-FM --output work/RNA01-02/rna-fm/CALIBRATION_REPEAT_1.npz
python paper_cases/rna/RNA01-02/models/rna-fm/forward/executed_full/encode.py --fasta work/RNA01-02/forward_input/CALIBRATION_UNIQUE_SEQUENCES.fasta --weights /path/to/rnafm_weights.pt --source-root /path/to/RNA-FM --output work/RNA01-02/rna-fm/CALIBRATION_REPEAT_2.npz
```

Run the public scorer using [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
