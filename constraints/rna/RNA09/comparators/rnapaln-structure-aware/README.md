# RNA09 RNApaln comparator

RNApaln is evaluated on the released RNA09 panel. The required
`--printAlignment` behavior is checked before the full run.

```bash
python constraints/rna/RNA09/comparators/rnapaln-structure-aware/run_rnapaln_primary.py \
  --panel constraints/rna/RNA09/refset \
  --rnapaln /path/to/RNApaln \
  --output-dir work/RNA09/rnapaln_primary

python constraints/rna/RNA09/comparators/rnapaln-structure-aware/score_rnapaln_primary.py \
  --panel constraints/rna/RNA09/refset \
  --score-dir work/RNA09/rnapaln_primary
```

## Released reference result

`results/RNAPALN_PRIMARY_RESULT.json` is the canonical comparator summary.
The corresponding raw directional measurements and ID-level F scores are
`results/RAW_DIRECTIONAL_SCORES.json` and
`results/RNAPALN_PRIMARY_LEAF_LEDGER.tsv`. Frozen parameters and the upstream
alignment regression check are stored beside them. The scoring program delegates
F-credit, hierarchy and bootstrap calculations to the RNA09 public scorer core.
