# RNA01-02 formal RefSet

This directory contains the complete, minimal, ready-to-Forward L16K release.

## Files

| file | role |
|---|---|
| `RNA0102_L16K_EDIT_REFSET/OBJECT_LEDGER.json` | 49,152 formal RNA objects and 49,149 unique sequence hashes with exact A/C/G/U sequences, role and grouping metadata; direct model Forward input |
| `RNA0102_L16K_EDIT_REFSET/EDIT_LADDER_LEDGER.json` | 16,384 frozen x0/x1/x2 ladders, object/hash mappings and outer-cluster identity |
| `RNA0102_L16K_EDIT_REFSET/PLANNED_POPULATION_MANIFEST.json` | frozen L16K membership and identity |
| `RNA0102_calibration_ladders.json` | 128 complete ladders used for two repeat Forward passes and scorer-internal `tau_delta` calibration |
| `CURRENT_RELEASE_STATUS.json` | public release status |

## Forward use

Encode every sequence in `OBJECT_LEDGER.json` once for formal scoring. Separately encode every endpoint referenced by `RNA0102_calibration_ladders.json` twice through the identical model path. Preserve `sequence_sha256` as the output row ID. The NPZ schema and scoring command are in [`../README.md`](../README.md).

No RefSet regeneration, candidate filtering or materialization pipeline is required.
