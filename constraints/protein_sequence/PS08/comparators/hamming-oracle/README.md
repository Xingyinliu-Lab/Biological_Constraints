# PS08 comparator: exact Hamming oracle

This deterministic truth-based oracle uses the same frozen PS08 RefSet and the same public formal scorer as foundation-model representations. It is a constructive positive control for the edit-ladder relation, not an independently trained traditional comparator and not a separately maintained scoring method.

## Method

For every frozen ladder `x0 -> x1 -> x2`, the program recomputes the four sequence Hamming distances and requires:

```text
d(x0,x1) = 1
d(x0,x2) = 2
d(x2,x1) = 1
d(x2,x0) = 2
```

It writes their negatives as similarity-like leaf scores and passes the resulting ledger to the single authoritative scorer:

```text
constraints/protein_sequence/PS08/scoring/score_ps08_formal.py
```

The comparator directory does not contain a second copy of that scorer.

## Run the comparator

From the repository root:

```bash
python constraints/protein_sequence/PS08/comparators/hamming-oracle/formalize_ps08.py \
  --oracle \
  --fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --outdir work/PS08/hamming-oracle
```

`--outdir` must not already exist. The script records the actual endpoint FASTA SHA256 and whether it matches the registered release file. It does not use byte identity as a general admission rule. It verifies the submitted population, sequence identities and every Hamming relation, creates `FORMAL_LEAF_SCORES.tsv`, invokes the public PS08 scorer with frozen `delta_num=0`, and writes `FORMAL_RESULT.json`.

## Outputs

| File | Contents |
|---|---|
| `FORMAL_LEAF_SCORES.tsv` | Frozen IDs and four Hamming-derived relation scores per ladder |
| `scorer/INDIVIDUAL_SCORES.tsv` | All 4,096 leaf IDs, raw measurements, F credits, hierarchy F scores and contributions |
| `scorer/SCORE_SUMMARY.json` | F, F 95% CI, aggregation counts and method identity |
| `scorer/BOOTSTRAP_F.float64le.bin` | Generated 20,000-draw F bootstrap array |
| `FORMAL_RESULT.json` | Comparator result plus FASTA and leaf-ledger identities |

The canonical flattened copies of these files, including
`INDIVIDUAL_SCORES.tsv`, are archived under [`results/`](results/) with an
updated `OUTPUT_MANIFEST.json`.
