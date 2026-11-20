# RNA04 released RefSet

This directory contains the minimal, ready-to-Forward RNA04 evaluation population.

| File | Role |
|---|---|
| `RNA04_unique_endpoints.fasta` | 11,935 unique coding-oriented WT/ALT sequences; FASTA IDs are sequence SHA256 values |
| `COMPLETE_SITE_ALT_LEDGER.tsv` | 11,506 formal variant rows with site, class, sequence, parent and dependency-cluster mappings |
| `COMPLETE_SITE_UNIVERSE_SUMMARY.json` | frozen population counts and materialization summary |
| `DEPENDENCY_CLOSURE.json` | dependency closure supporting the 554 formal outer clusters |

Encode every FASTA record and preserve its ID. Do not regenerate, filter or select a subset of this population. The Forward-output schema and scorer command are in the [RNA04 README](../README.md).

Upstream candidate manifests are outside the public release and are not required to evaluate a model.
