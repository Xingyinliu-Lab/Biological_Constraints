# DNA03 formal RefSet

This directory contains the complete, minimal and ready-to-Forward DNA03 evaluation population.

| File | Role |
|---|---|
| `DNA03_FORWARD_VIEWS.fasta` | 1,400 formal DNA view objects to encode |
| `FORWARD_VIEW_MANIFEST.json` | FASTA identity and object-count record |
| `CORE_VIEW_MEMBERSHIP.tsv` | view-to-core membership, dependency cluster and frozen core coordinates |
| `SEQUENCE_OBJECT.tsv` | formal object identity and sequence hashes |
| `PRIMARY_QUERY.tsv` | 1,120 frozen query/true-mate relations |
| `PRIMARY_NEGATIVE.tsv` | 35,840 frozen query/negative relations, 32 per query |
| `DNA03_TAU_CONTROL_IDS_v0.1.1.tsv` | 128 frozen repeat controls used to compute tau |
| `DNA03_REFSET_CONTRACT_v0.2.5_FROZEN_20261020.json` | current frozen RefSet contract |

Use `DNA03_FORWARD_VIEWS.fasta` directly; no candidate filtering or RefSet construction is required. Preserve all IDs and core coordinates. After Forward, save the representations according to the contract in [`../README.md`](../README.md) and pass them to the public scorer.

DNA03 v0.2.6 is the active F-only release. The released RefSet contract is [`DNA03_REFSET_CONTRACT_v0.2.5_FROZEN_20261020.json`](DNA03_REFSET_CONTRACT_v0.2.5_FROZEN_20261020.json).
