# PS01 formal RefSet

This directory is the ready-to-use formal input population for [PS01](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete promoted population. Do not substitute a historical small/S/4K set when a final large/L/16K population exists.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `CONFORMANCE_REPORT.json` | 740 | JSON object; top-level: `status`, `release_id`, `counts`, `validation` |
| `FORMAL_CELL_ALLOCATION_LEDGER.tsv.gz` | 52985 | columns: `role`, `clan_id`, `anchor_family`, `anchor_occurrence`, `positive_occurrence`, `negative_family`, `negative_occurrence`, `negative_parent`, `source_role_interpretation`, `allocation_status`, `candidate_negative_endpoints` |
| `FORMAL_CONTRAST_LEDGER.tsv.gz` | 86348 | columns: `role`, `clan_id`, `anchor_family`, `a`, `p`, `n`, `negative_family`, `anchor_parent`, `positive_parent`, `negative_parent`, `identity_p_num`, `identity_p_den`, `identity_n_num`, `identity_n_den`, `length_ratio_pn`, `composition_l1_5_pn`, `same_superkingdom_pn`, `superkingdom_anchor`, `superkingdom_positive`, `superkingdom_negative` |
| `FORMAL_OBJECT_LEDGER.tsv.gz` | 253853 | columns: `occurrence_id`, `role`, `clan_id`, `pfam_accession`, `parent_protein_id`, `stable_accession`, `sequence_version`, `superkingdom`, `species_taxon_id`, `sequence_sha256`, `source_sequence_sha256` |
| `FORMAL_SEQUENCES.faa.gz` | 187932 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `RELEASE_MANIFEST.json` | 2438 | JSON object; top-level: `release_id`, `files` |
| `SOURCE_SEAL.json` | 352 | JSON object; top-level: `status`, `source_candidate`, `historical_membership`, `membership_changed`, `formal_relation_build_authorized`, `formal_assay_refset_release` |

## Handoff to model execution

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline. Paper-specific checkpoint and forward provenance is maintained separately and does not alter the score calculation.
