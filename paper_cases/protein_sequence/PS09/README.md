# PS09 paper implementation cases

This directory records the five foundation-model implementations reported for PS09. Each case provides a configurable Forward entry and calls the single public scorer in `constraints/protein_sequence/PS09/scoring/score_ps09.py`.

The four non-ESM3 cases share the configurable Forward entry at [`shared_forward/ps09_four_model_forward.py`](shared_forward/ps09_four_model_forward.py). ESM3 retains its case-specific configurable residue extractor.

| Model ID | Case |
|---|---|
| `esm2-650m` | [`models/esm2-650m/`](models/esm2-650m/) |
| `prott5` | [`models/prott5/`](models/prott5/) |
| `esmc-600m` | [`models/esmc-600m/`](models/esmc-600m/) |
| `carp-640m` | [`models/carp-640m/`](models/carp-640m/) |
| `esm3-seq` | [`models/esm3-seq/`](models/esm3-seq/) |

The large endpoint-representation tensors are not published in this directory. Each Forward implementation regenerates them from the released PS09 RefSet and recorded checkpoint. Case outputs contain the public scorer's ID-level and aggregate F results.

PS09-specific scoring, cross-fit, folds and bootstrap rules are not duplicated here. They belong to the Constraint release.
