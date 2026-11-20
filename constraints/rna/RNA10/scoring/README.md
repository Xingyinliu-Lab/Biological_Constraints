# RNA10 scoring

`build_measurement_from_vectors.py` builds the public held-out measurement ledger from an NPZ and the released RefSet. `score_rna10.py` converts that ledger into F-only individual, block, source-unit and final outputs. Both programs default to `RNA10_FORMAL_RUNPLAN.json` and `RNA10_FORMAL_SCORING_BINDING.json` in this directory.

The formal outputs are `MEASUREMENT_LEDGER.json`, `INDIVIDUAL_SCORES.json` with raw decision scores and `F_credit`, `AGGREGATION_SCORES.json` with block/source F values, `FORMAL_RESULT.json` with `F` and `F_95_CI`, and `OUTPUT_MANIFEST.json`.

The measured model representations must exactly cover the RefSet sequence universe. All cross-fit, aggregation and bootstrap parameters are fixed by the supplied public scoring files.
