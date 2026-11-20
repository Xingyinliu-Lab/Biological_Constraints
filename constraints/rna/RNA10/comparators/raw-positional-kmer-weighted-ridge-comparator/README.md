# RNA10 positional k-mer ridge comparator

This deterministic sequence-only comparator creates positional k-mer vectors, builds a held-out RNA10 measurement ledger, and uses the shared F-only scorer.

From the repository root:

```bash
python constraints/rna/RNA10/comparators/raw-positional-kmer-weighted-ridge-comparator/code/build_features.py \
  constraints/rna/RNA10/refset/RNA10_UNIQUE_INPUTS.fasta \
  work/RNA10/comparator/representations.npz

python constraints/rna/RNA10/comparators/raw-positional-kmer-weighted-ridge-comparator/code/build_measurement.py \
  --effects constraints/rna/RNA10/comparators/raw-positional-kmer-weighted-ridge-comparator/method/EFFECT_EXECUTION_LEDGER.tsv \
  --vectors work/RNA10/comparator/representations.npz \
  --run-plan constraints/rna/RNA10/comparators/raw-positional-kmer-weighted-ridge-comparator/RNA10_COMPARATOR_RUNPLAN.json \
  --binding constraints/rna/RNA10/scoring/RNA10_FORMAL_SCORING_BINDING.json \
  --output work/RNA10/comparator/MEASUREMENT_LEDGER.json \
  --model-id raw_positional_kmer_weighted_ridge

python constraints/rna/RNA10/scoring/score_rna10.py \
  --run-plan constraints/rna/RNA10/comparators/raw-positional-kmer-weighted-ridge-comparator/RNA10_COMPARATOR_RUNPLAN.json \
  --measurement-ledger work/RNA10/comparator/MEASUREMENT_LEDGER.json \
  --model-id raw_positional_kmer_weighted_ridge \
  --out work/RNA10/comparator/FORMAL_RESULT.json
```

The released `results/` directory contains raw measurements and F-only individual, aggregate and formal score outputs.
