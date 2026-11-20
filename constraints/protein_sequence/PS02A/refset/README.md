# PS02A formal RefSet

This directory is the ready-to-use formal input population for [PS02A](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete promoted population. Do not substitute a historical small/S/4K set when a final large/L/16K population exists.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `CURRENT_RELEASE_STATUS.json` | 360 | JSON object; top-level: `constraint_id`, `status`, `supersedes_embedded_materialization_status`, `note` |
| `endpoint_sequences.faa` | 2051001 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `FILE_POINTERS.tsv` | 490 | columns: `role`, `path`, `sha256` |
| `FORMAL_BUNDLE_LEDGER.tsv.gz` | 1268553 | columns: `bundle_id`, `anchor`, `positive`, `negative`, `anchor_clan`, `anchor_family`, `anchor_parent`, `positive_clan`, `positive_family`, `positive_parent`, `negative_clan`, `negative_family`, `negative_parent`, `anchor_sequence_sha256`, `positive_match_M`, `positive_match_D`, `negative_match_M`, `negative_match_D` |
| `FORMAL_SOURCE_VALIDATION.json` | 368 | JSON object; top-level: `source_membership_status`, `formal_status`, `selection`, `matching`, `uncertainty`, `model_output_consumed` |
| `CURRENT_PUBLIC_MANIFEST.json` | 906 | Current consumable RefSet assets and identities |
| `POINT_DEPENDENCY_MEMBERSHIP.tsv.gz` | 48966 | columns: `bundle_id`, `group_id` |
| `PS02A_MATCHING_POLICY_AMENDMENT_v0.1.6_FROZEN_20261006.md` | 989 | human-readable frozen amendment or release note |
| `PS02A_UNCERTAINTY_AMENDMENT_v0.1.1_FROZEN_20261006.md` | 1363 | human-readable frozen amendment or release note |
| `RELEASE_VALIDATION.json` | 514 | JSON object; top-level: `release_id`, `formal_status`, `formal_refset_release`, `source_membership_validation`, `formal_relation_or_bundle_count`, `anchor_clan_bootstrap_units`, `point_membership`, `withdrawn_rule_artifact`, `model_output_consumed`, `model_scoring_status` |
| `UNCERTAINTY_ANCHOR_CLAN_MEMBERSHIP.tsv.gz` | 40599 | columns: `bundle_id`, `anchor_clan_bootstrap_unit` |

## Handoff to model execution

Encode all 13,720 records in `endpoint_sequences.faa` and preserve each FASTA
ID. The formal bundle ledger references 13,573 of those objects. The complete
13,720-object FASTA is nevertheless the frozen Forward input identity used by
the registered paper representations; users must not infer or materialize a
different subset. The public scorer validates exact FASTA coverage before using
only bundle-referenced objects in measurement calculations.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline.
