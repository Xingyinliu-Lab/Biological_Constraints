# DNA02-02 Hamming Level-0 reference

This frozen deterministic white-box reference reads the released sequences and verifies the exact edit truth

\[
H(x_0,x_1)=1,\quad H(x_1,x_2)=1,\quad H(x_0,x_2)=2.
\]

It is a geometric truth ceiling, not a learned competing model. Its canonical result is F=1, with 95% CI [1,1]. No additional traditional comparator is required.

Run from the repository root:

```bash
python constraints/dna/DNA02-02/comparators/hamming-level0/score_dna0202_hamming_level0.py \
  --output work/DNA02-02/hamming-level0
```

The program automatically loads the released FASTA and ladders, validates their frozen identities, and writes ID-level ladder measurements, query/locus/cluster ledgers, `FORMAL_RESULT.json` and `OUTPUT_MANIFEST.json`. Its ladder-to-query-to-locus-to-cluster aggregation and bootstrap are imported from the public scoring core. Canonical released outputs are retained in `results/` within this comparator directory.
