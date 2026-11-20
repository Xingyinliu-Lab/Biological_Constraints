# DNA06 formal RefSet

This directory is the complete ready-to-Forward DNA06 evaluation population.

| File | Role |
|---|---|
| `records.fasta` | 11,935 formal WT/ALT sequence objects |
| `sequence_object.tsv` | object IDs, sequence identities and lengths |
| `parent_cds.tsv` | 560 parent CDS objects and dependency clusters |
| `codon_site.tsv` | 3,920 formal codon sites and frozen coordinates |
| `variant_endpoint.tsv` | 11,506 single-nucleotide ALT endpoints and truth labels |
| `site_variant_membership.tsv` | exact site-to-variant relation membership |
| `FORWARD_TARGETS.tsv` | direct 15,426-row map from scorer row ID to released object and canonical codon interval |

No candidate filtering or RefSet construction is required. Use `FORWARD_TARGETS.tsv` to identify the sequence object and the 0-based half-open codon interval for each required codon-local row. Map model-specific tokens to those frozen biological coordinates during Forward, then submit the NPZ described in the Constraint README.

The released formal interface is `FORWARD_TARGETS.tsv` plus `../scoring/INPUT_CONTRACT.json`.
