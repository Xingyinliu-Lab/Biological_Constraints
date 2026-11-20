# PS05 comparator: NW plus BLOSUM62

This comparator uses global Needleman-Wunsch alignment with the fixed BLOSUM62 settings to produce the standard PS05 leaf ledger. Final F computation uses the PS05 public scorer.

## Produce the leaf ledger

```bash
python constraints/protein_sequence/PS05/comparators/nw-blosum62/produce_ps05_nw_blosum62.py \
  --relations constraints/protein_sequence/PS05/refset/FORMAL_RELATION_LEDGER.tsv \
  --plan constraints/protein_sequence/PS05/refset/PLANNED_SCOREABLE_POPULATION.json \
  --fasta constraints/protein_sequence/PS05/refset/PS05_ENDPOINT_SEQUENCES.faa \
  --tool <NEEDLE_EXECUTABLE> \
  --workers 32 \
  --outdir work/PS05/nw-blosum62/measurement
```

The output directory contains `FORMAL_LEAF_SCORES.tsv`, the pair-score cache and the comparator summary. All planned identities must be present.

## Score

```bash
python constraints/protein_sequence/PS05/scoring/score_ps05_formal.py \
  --contract-root constraints/protein_sequence/PS05/refset \
  --planned-population constraints/protein_sequence/PS05/refset/PLANNED_SCOREABLE_POPULATION.json \
  --leaf-ledger work/PS05/nw-blosum62/measurement/FORMAL_LEAF_SCORES.tsv \
  --delta-num 0 \
  --reuse-ledger constraints/protein_sequence/PS05/refset/ENDPOINT_REUSE_LEDGER.tsv \
  --individual-scores work/PS05/nw-blosum62/score/INDIVIDUAL_SCORES.tsv.gz \
  --output-dir work/PS05/nw-blosum62/score
```

The canonical published output is [`results/SCORE_SUMMARY.json`](results/SCORE_SUMMARY.json), which contains F, F 95% CI, coverage and method identity. ID-level measurements and aggregation outputs are regenerated locally from the released RefSet and comparator command.
