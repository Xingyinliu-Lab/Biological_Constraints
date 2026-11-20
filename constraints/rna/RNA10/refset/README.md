# RNA10 formal RefSet

This directory is the ready-to-use formal input population for [RNA10](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete RNA10 population: all unique WT and mutant sequences required by the 60 released matched blocks and 391 formal relations.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `block_members.jsonl` | 66861 | JSONL asset |
| `blocks.jsonl` | 36965 | JSONL asset |
| `C_RNA10_ASSAY_REFSET_CONTRACT_SCHEMA_v1.0.4_20261009.json` | 13624 | JSON object; top-level: `$schema`, `$id`, `type`, `additionalProperties`, `$defs`, `required`, `properties` |
| `MANIFEST.json` | 4908 | JSON object; top-level: `artifact`, `coverage`, `payloads`, `version` |
| `objects_mutant.jsonl` | 216321 | JSONL asset |
| `objects_wt.jsonl` | 25317 | JSONL asset |
| `RNA10_ASSAY_REFSET_CONTRACT_v1.0.4_FROZEN_20261009.json` | 19394 | JSON object; top-level: `artifact_type`, `schema_version`, `status`, `publication_status`, `constraint_id`, `parent_bindings`, `truth_semantics`, `reference_population`, `source_dependency_ledger`, `source_qc_gate`, `lifecycle_bindings`, `comparator_boundary`, `deferred_to_numeric_method`, `release_gate`, `prohibitions` |
| `RNA10_ASSAY_REFSET_RELEASE_v1.0.2_RETROACTIVE_20261009.json` | 13036 | JSON object; top-level: `artifact_type`, `release_id`, `status`, `publication_status`, `constraint_id`, `reference_design_release_ref`, `assay_refset_candidate_ref`, `reference_eligibility_payload_artifact_ref`, `reference_audit_report_release_ref`, `source_materialization_ref`, `population_identity_sha256`, `eligible_bundle_ids_hash`, `release_outputs`, `source_qc_evidence`, `gene_mapping_audit`, `candidate_semantic_content_hash`, `eligibility_ledger_hash` |
| `RNA10_UNIQUE_INPUTS.fasta` | 106887 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `source_dependency_units.jsonl` | 7380 | JSONL asset |
| `truth_relations.jsonl` | 153755 | JSONL asset |
| `uncertainty_dependency_clusters.jsonl` | 7860 | JSONL asset |
| `VALIDATION.json` | 631 | JSON object; top-level: `checks`, `status` |

## Handoff to model execution

Read every record in `RNA10_UNIQUE_INPUTS.fasta`, preserve its
`sequence_sha256` identifier and produce one whole-RNA vector per record.
`objects_wt.jsonl`, `objects_mutant.jsonl` and `truth_relations.jsonl` connect
those sequence identities to the released WT/mutant relations. Tokenization,
windowing and whole-sequence pooling belong to model Forward provenance.
RNA10 does not require a nucleotide-position representation.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline. Paper-specific checkpoint and forward provenance is maintained separately and does not alter the score calculation.
