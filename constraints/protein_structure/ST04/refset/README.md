# ST04 formal RefSet

This directory is the ready-to-use formal input population for [ST04](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete promoted population. Do not substitute a historical small/S/4K set when a final large/L/16K population exists.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `dependency.tsv` | 39234 | columns: `source_dependency_cluster_id`, `member_object_count`, `member_object_ids`, `dependency_basis` |
| `FORMAL_COORDINATE_INPUT_MANIFEST.json` | 80356 | JSON object; top-level: `object_count`, `records`, `status` |
| `manifest.json` | 2256 | JSON object; top-level: `artifact_id`, `constants`, `counts`, `files`, `source_processing_rule_id`, `status` |
| `objects/experimental_pdb/` | 50548954 | 586 released PDB coordinate assets; every file is referenced by one admitted formal object in `objects.tsv`. No additional source-materialization PDB is required for scoring. |
| `objects.tsv` | 182061 | columns: `structure_object_id`, `coordinate_locator`, `coordinate_sha256`, `canonical_sequence_sha256`, `canonical_sequence_length`, `coordinate_residue_count`, `coordinate_segment_count`, `parent_protein_id`, `source_dependency_cluster_id`, `mapping_method`, `mapping_cost`, `object_status` |
| `PUBLIC_REFSET_VALIDATION.json` | generated release record | Confirms exact relation-object closure and verifies all 586 coordinate hashes. |
| `residue_evidence.tsv` | 7012192 | Frozen mapping from each `residue_id` to `structure_object_id`, canonical polymer position and coordinate residue identity. It contains 98,209 data rows and has SHA256 `e0d1c4fa106959fa4ee1c142f0ec18576b806fcb8c7cf45f4a35df300e653c05`. |
| `relations.tsv` | 57544808 | columns: `relation_id`, `query_id`, `positive_edge_id`, `direction_id`, `structure_object_id`, `anchor_residue_id`, `true_target_residue_id`, `negative_residue_id`, `anchor_canonical_polymer_position`, `true_target_canonical_polymer_position`, `negative_canonical_polymer_position`, `anchor_insertion_code_normalization`, `true_target_insertion_code_normalization`, `negative_insertion_code_normalization`, `break_segment_eligibility`, `target_aa_identity`, `negative_aa_identity`, `aa_exact_match`, `signed_direction_class`, `signed_direction_source_positions`, `sequence_separation_bin`, `positive_distance`, `negative_distance`, `source_dependency_cluster_id`, `matching_status`, `within_query_rank`, `reuse_count`, `reuse_audit_status`, `execution documentation_id` |
| `ST04_ASSAY_REFSET_RELEASE_v0.1.1_FROZEN_20261012.json` | 1045 | JSON object; formal population counts and public scientific-table identities |
| `ST04_SPLIT_RELEASE_v0.1.1_FROZEN_20261012.json` | 118139 | JSON object; top-level: `algorithm`, `artifact_id`, `inner`, `inner_fold_counts_by_outer`, `inner_folds`, `outer`, `outer_fold_counts`, `outer_folds`, `relations_sha256`, `release_note`, `status`, `unit` |
| `VALIDATION.json` | 1155 | JSON object; top-level: `artifact_id`, `checks_executed`, `checks_failed`, `materialization_manifest_sha256`, `source_processing_rule_id`, `source_processing_rule_sha256`, `status`, `summary` |

## Handoff to model execution

Encode all 586 structures in `objects.tsv`. Use `coordinate_locator` as the
released coordinate input and preserve each `structure_object_id` when mapping
residue representations to the canonical polymer positions referenced by
`relations.tsv`. Use `residue_evidence.tsv` for the frozen coordinate-to-canonical
mapping rather than reconstructing that mapping from PDB numbering. The public package is already the exact object closure of the
152,180 formal relations; no source-materialization candidate needs to be
encoded.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline. Paper-specific checkpoint and forward execution documentation is maintained separately and does not alter the score calculation.
