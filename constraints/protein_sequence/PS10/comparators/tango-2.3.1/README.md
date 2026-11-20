# PS10 comparator: TANGO 2.3.1

TANGO 2.3.1 is the built-in sequence-based PS10 comparator. This directory contains the method, software environment, raw record scores, held-out ID-level results and canonical aggregate F result.

## Generate the frozen record scores

Follow [`ENVIRONMENT.md`](ENVIRONMENT.md), then run:

```bash
python constraints/protein_sequence/PS10/comparators/tango-2.3.1/run_tango_ps10.py \
  constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  /path/to/TANGO \
  work/PS10/tango/measurement
```

The producer uses the frozen pH, temperature, ionic-strength and terminal settings embedded in the released implementation. Their identity is recorded in [`method/TANGO_COMPARATOR_SPECIFICATION.json`](method/TANGO_COMPARATOR_SPECIFICATION.json).

## Apply the shared PS10 aggregation

```bash
python constraints/protein_sequence/PS10/comparators/tango-2.3.1/score_ps10_tango.py \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --ledger constraints/protein_sequence/PS10/refset/PRIMARY_ASSIGNMENTS.tsv \
  --scores work/PS10/tango/measurement/TANGO_SCORES.tsv \
  --numeric-method constraints/protein_sequence/PS10/scoring/PS10_NUMERIC_METHOD.json \
  --specification constraints/protein_sequence/PS10/comparators/tango-2.3.1/method/TANGO_COMPARATOR_SPECIFICATION.json \
  --output-dir work/PS10/tango/score
```

[`results/`](results/) contains the canonical TANGO record scores, held-out individual scores, fold scores, F bootstrap draws and `SCORE_SUMMARY.json`.

The producer validates the record schema, canonical peptide identity, class labels, lengths and sequence hashes. It records the actual input identities in its output.
