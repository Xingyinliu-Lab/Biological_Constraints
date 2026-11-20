# PS02B scoring

`score_ps02b.py` is the single public entry from model Forward output to formal
PS02B results.

Its internal stages are:

1. `representations_to_leaf.py`: exact Forward-FASTA coverage validation,
   ID alignment and frozen float32 total safe-cosine measurements;
2. `score_ps02b_formal.py`: relation credits, frozen hierarchy, dependency
   bootstrap, ID-level output and aggregate statistics.

These files are consecutive parts of one scorer, not alternative paper and
external-model implementations. Formal scientific parameters are constants and
cannot be changed through the public CLI.

See [`../README.md`](../README.md) for input schemas, commands, output columns
and interpretation. See [`ENVIRONMENT.md`](ENVIRONMENT.md) for dependencies.
