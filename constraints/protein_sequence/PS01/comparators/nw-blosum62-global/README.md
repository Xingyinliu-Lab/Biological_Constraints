# PS01 comparator: NW + BLOSUM62 global alignment

This is PS01's registered sequence-alignment reference. It uses the frozen PS01
RefSet and the same hierarchy, bootstrap unit, draw count and seed as the model
evaluation. It is not a separate constraint.

## Environment

Install EMBOSS `needle` and the Python dependencies recorded in
[`ENVIRONMENT.md`](ENVIRONMENT.md).

## Run the complete comparator

```bash
python constraints/protein_sequence/PS01/comparators/nw-blosum62-global/ps01_broad16k_nw_formal_wrapper_v2_0_1.py \
  --contrasts constraints/protein_sequence/PS01/refset/FORMAL_CONTRAST_LEDGER.tsv.gz \
  --fasta constraints/protein_sequence/PS01/refset/FORMAL_SEQUENCES.faa.gz \
  --needle /path/to/needle \
  --outdir work/PS01/nw_blosum62
```

`--workers` controls parallel execution and does not alter the formal population
or scoring definition. The wrapper produces the alignment-derived F-credit
ledger and the formal aggregate result. `ps01_nw_blosum62_comparator.py` is the
shared raw-alignment helper used by the wrapper.

The wrapper sends the completed contrast ledger to the same public F-only hierarchy and clan bootstrap used for foundation-model representations.

## Registered result

The canonical execution products are bundled with the Constraint:

- `results/RAW_MEASUREMENTS.tsv.gz`: alignment-derived pair measurements.
- `results/INDIVIDUAL_SCORES.tsv.gz`: the 1,256 formal F credits.
- `results/F_BOOTSTRAP.npy`: the 20,000 clan-bootstrap draws on the F scale.
- `results/FORMAL_RESULT.json`: aggregate F and its clan-bootstrap confidence interval.
- `results/COMPARATOR_LEDGER_SUMMARY.json`: population and ledger QC.

These outputs are reproducible from the released RefSet and comparator code;
they do not depend on `paper_cases/`.

