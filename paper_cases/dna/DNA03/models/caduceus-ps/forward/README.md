# DNA03 / caduceus-ps Forward

This directory preserves the model-specific code used to produce the archived DNA03 representations. It performs model Forward and core-coordinate extraction only; all scoring is delegated to the public Constraint scorer.

Formal input: `constraints/dna/DNA03/refset/DNA03_FORWARD_VIEWS.fasta` plus `CORE_VIEW_MEMBERSHIP.tsv`.

The published adapter emits `EMBEDDINGS.npz` and `REPEAT_EMBEDDINGS.npz` directly with the required `view_ids + embeddings` arrays. The core-coordinate extraction is unchanged from the executed implementation; generation of non-formal arrays is omitted.

The adapter records model-token spans (or a verified native one-nucleotide mapping), excludes special and padding tokens, averages all hidden states covering each nucleotide, and requires full nucleotide coverage. It uses the frozen 100/250/500/1000/2000-nt windows; for a 2,000-nt view it averages overlapping global nucleotide states from `[0,1000)`, `[500,1500)` and `[1000,2000)`, then mean-pools the published half-open core interval. Its Forward manifest and validation output are the mapping evidence; scorer input remains only `view_ids + embeddings`.

See [`../runplan/COMMANDS.md`](../runplan/COMMANDS.md) for the execution and scoring handoff.
