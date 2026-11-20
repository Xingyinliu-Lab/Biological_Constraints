# PS04 comparator

PS04 publishes NW + BLOSUM62 as its single formal comparator. It starts from
the released RefSet, constructs the raw sequence-alignment measurements, and
delegates relation aggregation and bootstrap to the unique PS04 public scorer.

| Comparator | F | 95% CI for F |
|---|---:|---|
| NW + BLOSUM62 | 0.237724 | [0.096656, 0.325520] |

The comparator directory contains executable code, environment instructions,
and canonical `results/`. Canonical
results include raw/leaf measurements, `INDIVIDUAL_SCORES.tsv.gz`,
`FORMAL_RESULT.json`, and `RESULT_MANIFEST.json`.

The EMBOSS-based method prescribes version 6.6.0.0 and fixed scoring
parameters, but records rather than requires the author's historical
executable SHA256.
