# ST01 frozen RefSet

This directory contains the complete, minimal released ST01 evaluation population. It is ready for model Forward and does not contain the upstream candidate structure pool or materialization intermediates.

## Model Forward assets

- `objects.tsv`: 28 formal object IDs, coordinate locators, coordinate hashes and dependency identities.
- `coordinates/`: the 28 exact PDB files to encode.

Encode all 28 objects. Do not replace coordinates, filter objects by model behavior or change the object IDs.
Coordinate filenames are locators, not formal object IDs: use `objects.tsv` to map each filename stem to its canonical `cath:...` ID before writing a global Forward NPZ.

## Scoring assets

- `relations.tsv`: 32 directed formal ranking contrasts.

## Built-in geometry reference

`pair_evidence.tsv` contains the frozen pairwise structural measurements required by the built-in Level-0 geometry reference. A third-party model does not read this file during Forward and must not use it to construct its representations.

After Forward, submit the NPZ described in the [ST01 README](../README.md) to the public scorer.
