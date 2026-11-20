# ST04 comparator: geometry reference

This deterministic Level-0 reference evaluates the released ST04 geometry
relations with the same F-only hierarchy as model scoring. It is a definitional
reference, not a trained predictor; its formal result is `F=1` and an empirical
confidence interval is not applicable.

From the repository root:

```bash
python constraints/protein_structure/ST04/comparators/geometry-whitebox/score_validate_st04_geometry_whitebox.py \
  --relations constraints/protein_structure/ST04/refset/relations.tsv \
  --result work/ST04/geometry-whitebox/FORMAL_RESULT.json \
  --validation work/ST04/geometry-whitebox/VALIDATION.json
```

To export the ID-level F measurements:

```bash
python constraints/protein_structure/ST04/comparators/geometry-whitebox/export_individual_scores.py \
  --relations constraints/protein_structure/ST04/refset/relations.tsv \
  --formal-result work/ST04/geometry-whitebox/FORMAL_RESULT.json \
  --output work/ST04/geometry-whitebox/INDIVIDUAL_SCORES.tsv
```
