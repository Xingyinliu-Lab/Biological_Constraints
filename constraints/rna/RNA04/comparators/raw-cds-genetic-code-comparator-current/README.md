# RNA04 raw-CDS genetic-code reference

This built-in Level-0 reference derives the single changed codon directly from complete WT/ALT CDS sequences using translation table 1. It assigns endpoint score 1 to an amino-acid-changing variant and 0 to a synonymous variant, then reuses the same public RNA04 site/parent/cluster aggregation core as model scoring.

It is a deterministic scientific reference owned by RNA04, not a foundation model and not a second scorer. Its pair, site, parent and cluster outputs use the same F-only aggregation scale as model scoring.

## Run

The released full-chain entry reads the bundled RefSet and writes endpoint, individual and aggregate outputs:

```bash
python constraints/rna/RNA04/comparators/raw-cds-genetic-code-comparator-current/run_comparator.py \
  --output work/RNA04/raw-cds-reference
```

`rna04_raw_cds_genetic_code_comparator.py` is the inspectable raw WT/ALT CDS endpoint method. `run_comparator.py` is the only complete public execution entry and delegates formal aggregation to `scoring/rna04_public_core.py`; no second aggregation implementation is maintained here.

Canonical reference outputs are in `results/`, including ID-level endpoint and site scores plus the formal aggregate result and identities.
