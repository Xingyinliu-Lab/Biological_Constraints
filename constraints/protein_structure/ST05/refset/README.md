# ST05 formal RefSet

This directory is the ready-to-use formal input population for [ST05](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete promoted population. Do not substitute a historical small/S/4K set when a final large/L/16K population exists.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `attrition.tsv` | 23 | columns: `object_id`, `level`, `reason` |
| `coordinates/` | 88659760 | 1020 released PDB coordinate assets; object-to-file mapping is defined by the RefSet manifest/table. |
| `dependency.tsv` | 2717 | columns: `base_object_id`, `source_parent_id`, `source_dependency_cluster_id` |
| `manifest.json` | 1405 | JSON object; top-level: `artifact_id`, `status`, `counts`, `payloads` |
| `objects.tsv` | 21293 | columns: `base_object_id`, `coordinate_locator`, `coordinate_sha256`, `coordinate_scope_sha256`, `coordinate_record_count`, `source_parent_id`, `source_dependency_cluster_id`, `canonical_sequence_sha256`, `canonical_sequence_length`, `smoke_topology_id`, `source_inventory_object_status`, `truth_status`, `provenance_id` |
| `relations.tsv` | 687579 | columns: `relation_id`, `panel_bundle_id`, `base_object_id`, `variant_object_id`, `transform_id`, `is_identity_transform`, `rotation_matrix`, `translation_vector`, `composition_order`, `generator_id`, `generator_version`, `seed`, `sampling_stratum`, `base_coordinate_sha256`, `variant_coordinate_sha256`, `coordinate_scope_sha256`, `source_parent_id`, `source_dependency_cluster_id`, `proper_rotation_status`, `uniform_translation_status`, `identity_order_scope_status`, `internal_distance_preservation_status`, `max_component_serialization_error_angstrom`, `all_pair_distance_error_bound_angstrom`, `variant_coordinate_locator`, `truth_status`, `provenance_id` |
| `ST05_ASSAY_REFSET_RELEASE_v0.1.0_FROZEN_20261012.json` | 766 | JSON object; formal population counts and public table identities |
| `transforms.tsv` | 4181 | columns: `transform_id`, `sampling_stratum`, `is_identity_transform`, `axis`, `angle_degrees`, `rotation_matrix`, `translation_vector`, `composition_order`, `generator_id`, `generator_version`, `seed` |

## Handoff to model execution

`objects.tsv` maps each `base_object_id` to its base PDB under
`coordinates/base/`; `relations.tsv` maps every `relation_id` and
`variant_object_id` to its transformed PDB under `coordinates/variants/`.
Models must encode the 60 bases and all 960 listed variants with the same
atom/residue selection, masking, preprocessing and global readout. The paper
cases provide a staging command that flattens this released panel while
preserving coordinate bytes and the required `<base>__BASE` / variant ID names.
Record those model-specific choices with the Forward output.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public F-only scoring pipeline.
