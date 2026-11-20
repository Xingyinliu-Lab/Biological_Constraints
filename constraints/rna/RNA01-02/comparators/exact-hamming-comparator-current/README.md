# RNA01-02 exact-Hamming reference

This is the frozen Level-0 white-box reference and deterministic ceiling for RNA01-02. It uses the defining exact edit relation

`H(x0,x1)=1`, `H(x1,x2)=1`, `H(x0,x2)=2`

with oracle `tau_delta=0`. It is an internal scientific reference for this Constraint, not a learned model or an independent traditional predictor.

## Reproduce the reference

From the repository root:

```bash
python constraints/rna/RNA01-02/comparators/exact-hamming-comparator-current/rna0102_hamming_summary_runner.py \
  --output constraints/rna/RNA01-02/comparators/exact-hamming-comparator-current/results
```

The runner automatically loads the released L16K ladder and object ledgers, validates sequence hashes and exact 1/1/2 Hamming relations, writes ID-level measurements, aggregates the same outer-cluster units, and verifies the deterministic ceiling.

## Outputs

- `results/INDIVIDUAL_SCORES.tsv`: one row per L16K ladder with object IDs, exact distances and endpoint F credits.
- `results/CLUSTER_SCORE_LEDGER.tsv`: dependency-closed `cluster_F` scores.
- `results/FORMAL_RESULT.json`: `F=1` and `F_95_CI=[1,1]`.
- `results/OUTPUT_MANIFEST.json`: output SHA256 and sizes.

`rna0102_hamming_comparator.py` contains the exact method and input validation. The summary runner is the public full-chain entry. No file under `paper_cases/` is required.
