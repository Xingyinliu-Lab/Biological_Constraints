# ST01 deterministic geometry reference

This built-in Level-0 white-box reference reads the frozen TM-align pair evidence, validates the 28 released coordinate files and applies the same ST01 relation and aggregation hierarchy. It is a deterministic ceiling/reference, not a learned structure predictor.

## Reproduce the reference

```bash
python constraints/protein_structure/ST01/comparators/geometry-comparator/score_validate_st01_geometry_level0.py \
  --relations constraints/protein_structure/ST01/refset/relations.tsv \
  --pair-evidence constraints/protein_structure/ST01/refset/pair_evidence.tsv \
  --objects constraints/protein_structure/ST01/refset/objects.tsv \
  --coordinate-root constraints/protein_structure/ST01/refset/coordinates \
  --numeric-method constraints/protein_structure/ST01/scoring/ST01_NUMERIC_METHOD_v0.1.1_FROZEN_20261010.json \
  --binding constraints/protein_structure/ST01/scoring/ST01_NUMERIC_INPUT_BINDING.json \
  --result work/ST01/geometry/REPRODUCED_RESULT.json \
  --validation work/ST01/geometry/VALIDATION.json \
  --individual-scores work/ST01/geometry/INDIVIDUAL_SCORES.tsv
```

The comparator is implemented in [`score_validate_st01_geometry_level0.py`](score_validate_st01_geometry_level0.py). Canonical public outputs are in [`results/`](results/): the formal result is `F=1`; uncertainty is not applicable to this definitional deterministic reference. `INDIVIDUAL_SCORES.tsv` exposes all 32 directed contrast margins and F credits.
