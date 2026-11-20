# PS02A scoring

`score_ps02a.py` is the only public user entry. It accepts model Forward output
and writes both ID-level measurements/scores and aggregate statistics.

Internal stages:

1. `representations_to_leaf.py` validates and aligns representations, then
   computes frozen float32 total safe-cosine measurements for every bundle.
2. `score_ps02a_formal.py` converts those measurements to formal relation
   credits, performs the frozen hierarchy and anchor-clan bootstrap, and writes
   the aggregate result.

The two internal files are consecutive stages of one public scorer, not
alternative scoring interfaces. Paper models and external models call
`score_ps02a.py` identically. Scientific parameters are constants in the formal
implementation and are not user options.

See [`../README.md`](../README.md) for the complete representation contract,
commands, output columns and interpretation. See [`ENVIRONMENT.md`](ENVIRONMENT.md)
for dependencies.
