# PS04 Needleman-Wunsch BLOSUM62 comparator

`ps04_direct_alignment_comparator.py` produces the pair-score assets and
`run_ps04_direct_alignment_formal.py` builds the formal comparator result. Final
population validation, aggregation and bootstrap use the unique PS04
implementation under [`../../scoring/`](../../scoring/README.md).

The comparator uses EMBOSS Needle 6.6.0.0 and the public PS04 F-only scorer.
Run it from the repository root with explicit local input and output paths:

```bash
python constraints/protein_sequence/PS04/comparators/nw-blosum62/run_ps04_direct_alignment_formal.py \\
  --tool <NEEDLE_EXECUTABLE> \\
  --outdir work/PS04/nw-blosum62
```

The runner fixes global Needleman-Wunsch with EBLOSUM62 and reads its required
RefSet and formal-scoring inputs from this public Constraint package. The
committed result records formal `F` and `F_95_CI`. Full pair, leaf and
relation-level ledgers are generated into the caller's output directory and
are not committed to the repository.
