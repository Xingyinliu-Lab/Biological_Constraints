# RNA01-02 × RiNALMo-Giga Forward

This directory preserves the Forward implementation used to produce the registered RNA01-02 representations. It does not contain a scoring algorithm.

## Inputs and outputs

The Forward input and representation contract are defined by the [RNA01-02 Constraint](../../../../../../constraints/rna/RNA01-02/README.md). The implementation must export the full L16K vectors and two calibration-repeat files with the exact object IDs required by that contract.

## Implementation

- Executed FULL-package entry: [`executed_full/encode.py`](executed_full/encode.py)
- Recovered supporting scripts: [`scripts/`](scripts/)
- Model identity and environment: [`paper_cases/model_resources/rna/rinalmo-giga/`](../../../../../model_resources/rna/rinalmo-giga/)
- Public scorer invocation: [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md)

The public Constraint scorer performs calibration, relation measurement, aggregation and confidence-interval calculation. Result-local scoring implementations are not part of the release.

## Current runnable sequence

From the repository root, first run the RefSet materialization command in the Constraint README. Then use the same extraction settings for the L16K population and for two independent calibration passes:

```bash
mkdir -p work/RNA01-02/rinalmo-giga
python paper_cases/rna/RNA01-02/models/rinalmo-giga/forward/executed_full/encode.py --fasta work/RNA01-02/forward_input/L16K_UNIQUE_SEQUENCES.fasta --repo /path/to/RiNALMo --out work/RNA01-02/rinalmo-giga/L16K_VECTORS.npz
python paper_cases/rna/RNA01-02/models/rinalmo-giga/forward/executed_full/encode.py --fasta work/RNA01-02/forward_input/CALIBRATION_UNIQUE_SEQUENCES.fasta --repo /path/to/RiNALMo --out work/RNA01-02/rinalmo-giga/CALIBRATION_REPEAT_1.npz
python paper_cases/rna/RNA01-02/models/rinalmo-giga/forward/executed_full/encode.py --fasta work/RNA01-02/forward_input/CALIBRATION_UNIQUE_SEQUENCES.fasta --repo /path/to/RiNALMo --out work/RNA01-02/rinalmo-giga/CALIBRATION_REPEAT_2.npz
```

Run the public scorer using [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md).
