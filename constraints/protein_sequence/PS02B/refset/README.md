# PS02B formal RefSet

This directory is the ready-to-use formal input population for [PS02B](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete promoted population. Do not substitute a historical small/S/4K set when a final large/L/16K population exists.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `CURRENT_RELEASE_STATUS.json` | 336 | JSON object; top-level: `constraint_id`, `status`, `supersedes_embedded_materialization_status`, `note` |
| `endpoint_sequences.faa` | 682411 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `FILE_POINTERS.tsv` | 615 | columns: `role`, `path_or_identity`, `sha256` |
| `FORMAL_OBJECT_LEDGER.tsv.gz` | 230946 | columns: `object_id`, `domain_id`, `cath_release`, `topology_id`, `hsf_id`, `parent_proxy_pdb_chain`, `sequence_length`, `sequence_sha256`, `domain_length_cath`, `resolution_or_quality` |
| `FORMAL_RELATION_LEDGER.tsv.gz` | 328217 | columns: `bundle_id`, `anchor`, `positive`, `negative`, `topology`, `anchor_hsf`, `negative_hsf`, `negative_parent`, `identity_pos`, `identity_neg`, `coverage_anchor_pos`, `coverage_anchor_neg`, `coverage_partner_pos`, `coverage_partner_neg` |
| `FORMAL_SCORING_PRECHECK.json` | 375 | JSON object; top-level: `membership_reconstructed_identically`, `source_rows`, `formal_rows`, `model_forward_required`, `formal_scoring_status`, `reason` |
| `FORMAL_SOURCE_VALIDATION.json` | 674 | JSON object; top-level: `release_id`, `formal_status`, `source_selected_relation_ledger`, `validation`, `note` |
| `MANIFEST.json` | 1693 | JSON object; top-level: `release_kind`, `formal_status`, `files` |
| `POINT_DEPENDENCY_MEMBERSHIP.tsv.gz` | 23036 | columns: `bundle_id`, `group_id` |
| `RELEASE_VALIDATION.json` | 571 | JSON object; top-level: `release`, `manifest_entries_before_finalization`, `manifest_errors_before_finalization`, `formal_ledger`, `formal_relation_or_bundle_count`, `group_counts`, `formal_status`, `scoring_status`, `note` |
| `UNCERTAINTY_GROUP_MEMBERSHIP.tsv.gz` | 23037 | columns: `bundle_id`, `group_id` |

## Handoff to model execution

Forward every record in `endpoint_sequences.faa`. The file contains 3,937
unique domain IDs. The formal relations score 3,905 of them; the remaining 32
are still part of the frozen released Forward population and must not be
silently omitted from a submitted representation set. Preserve each FASTA ID
exactly and emit one whole-domain vector per record as described in the main
PS02B README.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline.
