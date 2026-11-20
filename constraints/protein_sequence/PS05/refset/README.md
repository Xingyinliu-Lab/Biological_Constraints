# PS05 formal RefSet

This directory is the ready-to-use formal input population for [PS05](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode all 321 objects in `PS05_ASSAY_REFSET_RELEASE_v1.0.1_20261008`.
  Historical development populations are not public alternatives.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `AggregationGraph.json` | 483 | JSON object; top-level: `contract`, `schema_id`, `schema_version`, `status` |
| `ATCHLEY_2005_FACTORS.tsv` | 743 | columns: `amino_acid`, `factor1`, `factor2`, `factor3`, `factor4`, `factor5` |
| `ENDPOINT_REUSE_LEDGER.tsv` | 7725 | columns: `endpoint_id`, `anchor_count`, `positive_count`, `negative_count`, `reuse_count`, `high_reuse_gt_100` |
| `FORMAL_RELATION_LEDGER.tsv` | 333759 | columns: `anchor_domain_id`, `positive_domain_id`, `negative_domain_id`, `anchor_hsf_id`, `positive_hsf_id`, `negative_hsf_id`, `topology_id`, `anchor_parent_protein_id`, `positive_parent_protein_id`, `negative_parent_protein_id`, `S_struct_anchor_positive`, `S_struct_anchor_negative`, `seqid_anchor_positive`, `seqid_anchor_negative`, `anchor_coverage_positive`, `anchor_coverage_negative`, `partner_coverage_positive`, `partner_coverage_negative` |
| `MetricRegistry.json` | 514 | JSON object; top-level: `contract`, `schema_id`, `schema_version`, `status` |
| `PLANNED_SCOREABLE_POPULATION.json` | 338232 | JSON object; top-level: `kind`, `constraint_id`, `assay_id`, `measurement_contract_ref`, `identity_fields`, `leaf_identities`, `scoreable_population_identity_sha256` |
| `PS05_ENDPOINT_SEQUENCES.faa` | 46680 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `ResultContract.json` | 612 | JSON object; top-level: `contract`, `schema_id`, `schema_version`, `status` |

## Handoff to model execution

Forward every FASTA record exactly once and retain the FASTA identifier as the submitted `object_id`. PS05 uses whole-domain representations, so it does not prescribe residue-to-token indices. Users are responsible for documenting their tokenizer, special-token handling, truncation or windowing, and whole-domain pooling rule.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline. Paper-specific checkpoint and forward provenance is maintained separately and does not alter the score calculation.
