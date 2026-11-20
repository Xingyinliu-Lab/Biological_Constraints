# ST02 comparator: TM-align

TM-align is the registered structural-alignment comparator for ST02. It measures
the released pair panel, then applies the same public F-only hierarchy and
20,000-draw dependency-cluster bootstrap as the model scorer.

## Environment

Install a compatible TM-align executable and the Python dependencies listed in
[`ENVIRONMENT.md`](ENVIRONMENT.md).

## Run

From the repository root, generate the raw pair measurements:

```bash
python constraints/protein_structure/ST02/comparators/tmalign-comparator/run_tmalign_pair_panel.py \
  --tmalign /path/to/TMalign \
  --output work/ST02/tmalign_pairs.json
```

Then create the formal F-only result and ledgers:

```bash
python constraints/protein_structure/ST02/comparators/tmalign-comparator/score_st02_tmalign.py \
  --pair-ledger work/ST02/tmalign_pairs.json \
  --output work/ST02/tmalign/FORMAL_COMPARATOR_RESULT.json
```

The scorer writes `F`, `F_95_CI`, the individual `F_credit` ledger, hierarchy
F ledgers, and `RESULT_MANIFEST.json` next to the requested result JSON. The
raw pair ledger records the TM-align executable identity and contains the
reproducible symmetric bidirectional mean TM-score measurements.
