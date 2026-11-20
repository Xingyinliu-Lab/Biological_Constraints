# RNA05 ViennaRNA pairing-propensity comparator

This is the frozen RNA05 comparator on its explicitly compatible scope: 1,359 of the 1,514 formal blocks. Each compatible block contributes one contrast: its mean paired-position propensity versus its mean explicitly-unpaired-position propensity. The remaining 155 blocks are reported as out of scope rather than silently removed from the full RNA05 population.

## Released assets

- `vienna_partition_pairing_propensity.py` computes per-nucleotide pairing propensity from the released RNA sequences with ViennaRNA.
- `results/RNA05_VIENNA_SCORE_LEDGER.json` preserves the 1,359 registered block-level raw measurements and F credits.
- `RECONSTRUCT_CURRENT_VIENNA_RESULT.py` validates the compatibility audit and score ledger, then recomputes compatible-scope F.
- `results/FORMAL_RESULT.json` is the registered comparator result, including F and F 95% CI.
- `results/RNA05_VIENNA_COMPATIBILITY_AUDIT.json` records the compatible and out-of-scope populations.

The executed environment is recorded in [`ENVIRONMENT.md`](ENVIRONMENT.md). Comparator results are not copied into `paper_cases/`.

## Reconstruct the registered aggregate

```bash
mkdir -p work/RNA05/vienna
python constraints/rna/RNA05/comparators/viennarna-pairing-propensity/RECONSTRUCT_CURRENT_VIENNA_RESULT.py \
  --score-ledger constraints/rna/RNA05/comparators/viennarna-pairing-propensity/results/RNA05_VIENNA_SCORE_LEDGER.json \
  --compatibility-audit constraints/rna/RNA05/comparators/viennarna-pairing-propensity/results/RNA05_VIENNA_COMPATIBILITY_AUDIT.json \
  --result constraints/rna/RNA05/comparators/viennarna-pairing-propensity/results/FORMAL_RESULT.json \
  --out work/RNA05/vienna/reconstructed_result.json
```

## Re-run from released RNA sequences

The portable adapter implements the frozen block-level contrast and uses no
author-machine paths:

```bash
python constraints/rna/RNA05/comparators/viennarna-pairing-propensity/code/run_vienna_comparator.py \
  --output-dir work/RNA05/vienna/full_rebuild
```

The command reads the released RNA05 FASTA and endpoint ledger plus the frozen
compatibility pair map under `method/`. It writes the compatibility audit,
ID-level raw measurements and F credits, and `FORMAL_RESULT.json`.

The registered result is reproduced by rerunning this command with the released
released assets; verify the resulting compatibility audit, score ledger and
formal result before replacing a registered result.
