# RNA09 formal RefSet

This directory is the ready-to-use formal input population for [RNA09](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode every unique sequence in `scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta`; the object map links sequence hashes to the formal endpoint objects.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `ANCHOR_BUNDLE_LEDGER.tsv` | 2951 | columns: `anchor_id`, `source_parent_id`, `source_dependency_unit_id`, `contrast_ids`, `released_contrast_count`, `complete_bundle_status` |
| `FORMAL_CONTRAST_LEDGER.tsv` | 13432 | columns: `contrast_id`, `anchor_id`, `source_parent_id`, `source_dependency_unit_id`, `positive_partner_id`, `negative_partner_id`, `positive_relation_id`, `negative_relation_id`, `length_ratio_band`, `candidate_rank_sha256`, `selection_rank_global`, `positive_role_reuse_count`, `negative_role_reuse_count`, `allocation_contract_id` |
| `INDEPENDENT_VALIDATION_SUMMARY.json` | 807 | JSON object; top-level: `allocation_independently_recomputed`, `artifact_type`, `assay_induced_reuse_excluded_from_source_dependency`, `cap_respected`, `counts`, `exact_selected_contrast_membership_match`, `planned_population_identity_reconstructed`, `role_specific_cap`, `source_dependency_reconstructed_from`, `source_unit_membership_match`, `status` |
| `MANIFEST.json` | 1220 | JSON object; top-level: `package`, `payload` |
| `MATERIALIZATION_SUMMARY.json` | 706 | JSON object; top-level: `artifact_type`, `cap`, `counts`, `independent_validation`, `note`, `planned_population_identity_sha256`, `reuse`, `scientific_contract_version`, `status` |
| `PLANNED_POPULATION_MANIFEST.json` | 1362 | JSON object; top-level: `allocation_contract_id`, `artifact_type`, `authority_model`, `contract_id`, `formal_bundle`, `input_hashes`, `planned_counts`, `planned_population_identity_sha256`, `prohibited_dependency_edges`, `relation_truth`, `scientific_contract_version`, `selection_rule`, `status` |
| `REUSE_LEDGER.tsv` | 2938 | columns: `role`, `object_id`, `released_reuse_count`, `frozen_cap`, `allocation_contract_id` |
| `scoring_input_v01/RNA09_ENDPOINT_OBJECT_MAP.tsv` | 3137 | columns: `sequence_sha256`, `object_id`, `sequence_length_nt` |
| `scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta` | 5079 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `scoring_input_v01/SUMMARY.json` | 371 | JSON object; top-level: `artifact_type`, `endpoint_objects`, `fasta_bytes`, `fasta_header`, `fasta_sha256`, `object_map_sha256`, `status`, `unique_sequence_hashes` |
| `SOURCE_DEPENDENCY_LEDGER.tsv` | 5539 | columns: `source_dependency_unit_id`, `member_object_ids`, `edge_reason_types` |

## Handoff to model execution

Run the user's model on the complete FASTA and retain each FASTA
`sequence_sha256` as the representation row ID. Tokenization, windowing and
whole-RNA pooling must be documented with the model Forward implementation; RNA09 does not require a
position-level nucleotide mapping.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline. Paper-model checkpoint and Forward details do not alter the score calculation.
