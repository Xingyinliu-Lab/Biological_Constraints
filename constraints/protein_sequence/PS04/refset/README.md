# PS04 formal RefSet

This directory is the ready-to-use formal input population for [PS04](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode all 3,641 objects in `FORMAL_ENDPOINTS.faa`. The unique active PS04
  population is `PS04_ASSAY_REFSET_RELEASE_v1.0.1_L_4K_20261007`; PS04 has no
  separate L/16K public population.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `CANDIDATE_LINEAGE.json` | 1163 | JSON object; top-level: `artifact_id`, `candidate_manifest_sha256`, `candidate_path`, `candidate_payload_sha256`, `formalization` |
| `ENDPOINT_REUSE_LEDGER.tsv.gz` | 33779 | columns: `role`, `endpoint_id`, `anchor_bundle_reuse_count` |
| `FILE_POINTERS.tsv` | 616 | columns: `role`, `path`, `sha256` |
| `FORMAL_BUNDLE_LEDGER.tsv.gz` | 90082 | columns: `bundle_id`, `anchor_id`, `anchor_architecture_id`, `anchor_topology_id`, `anchor_hsf_id`, `anchor_parent_proxy`, `anchor_source_unit_id`, `legal_leaf_count`, `legal_positive_count`, `legal_negative_hsf_count`, `status` |
| `FORMAL_ENDPOINTS.faa` | 585483 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `FORMAL_OBJECT_LEDGER.tsv.gz` | 305142 | columns: `rank`, `object_id`, `domain_id`, `cath_release`, `architecture_id`, `topology_id`, `hsf_id`, `parent_proxy_pdb_chain`, `source_unit_id`, `sequence_sha256`, `sequence_length` |
| `FORMAL_RELATION_LEDGER.tsv.gz` | 60354295 | columns: `anchor_id`, `positive_id`, `negative_id`, `negative_hsf_id`, `negative_topology_id`, `anchor_architecture_id`, `anchor_topology_id`, `anchor_hsf_id`, `anchor_parent_proxy`, `anchor_source_unit_id`, `positive_hsf_id`, `positive_parent_proxy`, `positive_source_unit_id`, `negative_parent_proxy`, `negative_source_unit_id`, `positive_seqid`, `negative_seqid`, `positive_anchor_coverage`, `negative_anchor_coverage`, `positive_partner_coverage`, `negative_partner_coverage` |
| `FORMAL_SOURCE_VALIDATION.json` | 1010 | JSON object; top-level: `anchor_hsfs`, `architectures`, `artifact_id`, `candidate_source_panel`, `complete_bundles`, `endpoint_reuse_summary`, `formal_source_membership_validated`, `leaf_relations`, `object_count`, `policy`, `source_dependency_units`, `status`, `strict_validation_failures`, `topologies` |
| `INPUT_PINS.json` | 657 | JSON object; top-level: `method`, `model_derived_inputs`, `negative_pair_ledger_sha256`, `positive_pair_ledger_sha256`, `source_object_table_sha256`, `source_traversal_sha256` |
| `MANIFEST.json` | 2802 | JSON object; top-level: `artifact_id`, `payloads`, `status` |
| `POINT_DEPENDENCY_MEMBERSHIP.tsv.gz` | 58382 | columns: `bundle_id`, `anchor_id`, `point_dependency_unit_id`, `reason_codes` |
| `SOURCE_DEPENDENCY_MEMBERSHIP.tsv.gz` | 71213 | columns: `object_id`, `source_unit_id`, `dependency_reason_codes` |
| `UNCERTAINTY_ANCHOR_SOURCE_UNIT_MEMBERSHIP.tsv.gz` | 58398 | columns: `bundle_id`, `anchor_id`, `anchor_source_unit_bootstrap_unit`, `reason_codes` |

## Handoff to model execution

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline.
