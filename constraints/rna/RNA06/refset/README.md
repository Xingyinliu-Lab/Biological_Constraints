# RNA06 RefSet

This directory is the complete, minimal, ready-to-Forward primary RNA06 population.

| Asset | Role |
|---|---|
| `rna_sequences.fasta` | 105 model-input RNA sequences keyed by sequence SHA256 |
| `rna_objects.tsv` | record-to-sequence and coordinate identity mapping |
| `RNA06_FORMAL_DIRECTED_ENDPOINT_LEDGER.tsv` | 3,176 directed true-partner challenges |
| `RNA06_FORMAL_PAIRED_NONPARTNER_CANDIDATE_LEDGER.tsv` | 14,948 frozen paired-nonpartner candidates |
| `RNA06_FORMAL_RECIPROCAL_EDGE_BUNDLE_LEDGER.tsv` | 1,588 reciprocal edge bundles and hierarchy IDs |
| `RNA06_FROZEN_PAIR_STATE_MAPS.tsv` | frozen positive-pair identities used to validate true partners and paired nonpartners |
| `PLANNED_POPULATION_MANIFEST.json` | frozen population identity and counts |

The FASTA contains each unique RNA sequence once and uses its sequence SHA256
as the record ID. `rna_objects.tsv` maps the 114 biological RNA object records
to these 105 unique model-input sequences.

The Vienna-compatible subset is comparator-specific and is stored under `comparators/viennarna-partner-probability/refset/`; it must not be used as the formal model population.

All ledger positions are zero-based biological nucleotide coordinates. Model token indices are not interchangeable with these positions. Users must document token-to-nucleotide span mapping, windowing and position-vector readout, and must resolve anchor, true-partner and nonpartner targets from their explicit sequence/coordinate fields rather than row order. See the [RNA06 README](../README.md) for the complete Mapping Evidence requirement.
