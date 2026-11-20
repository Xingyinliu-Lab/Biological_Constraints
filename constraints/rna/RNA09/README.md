# RNA09 | Remote RNA Structure

RNA09 tests whether whole-RNA representations distinguish a structurally similar partner from a matched, lower-similarity alternative under the released remote-sequence and length-ratio design.

## Public contents

`refset/` contains the formal endpoint sequences, object map, contrast ledger, anchor bundles, dependency units, reuse ledger and population manifest. `scoring/RNA09_NUMERIC_METHOD_v0.1.3_F_ONLY.json` defines the formal F-only score. `comparators/` contains the RNApaln comparator and its released result.

## Representation input

Provide one whole-RNA vector per released endpoint sequence in an NPZ with exactly two arrays: `sequence_sha256` and `vectors`. The IDs must exactly match the released endpoint FASTA. Model-specific tokenization, special-token handling, windowing and pooling must produce one finite vector per sequence.

## Score a new model

From the repository root:

```bash
python constraints/rna/RNA09/scoring/score_official.py \
  --embeddings work/RNA09/<model>/REPRESENTATIONS.npz \
  --output-dir work/RNA09/<model>/score
```

The scorer validates the released RefSet identity and emits F-only leaf, anchor, parent, dependency-unit and group ledgers, `F_BOOTSTRAP_MEANS.npy`, `FORMAL_RESULT.json` and `OUTPUT_MANIFEST.json`.

## Formal score

Each contrast receives `F_credit = 1` when the positive similarity exceeds the matched negative similarity, `0` for an exact tie and `-1` otherwise. Equal immediate-child means aggregate through:

`negative -> positive -> anchor -> source parent -> source dependency unit -> structure relation group -> F`.

`F_95_CI` is the two-sided percentile interval from 20,000 dependency-unit bootstrap draws (seed 20261004). The score describes performance on this released reference population; it is not a universal measure of RNA structural understanding.

## Comparator

See [`comparators/rnapaln-structure-aware/README.md`](comparators/rnapaln-structure-aware/README.md) for RNApaln measurement and F-only scoring commands.
