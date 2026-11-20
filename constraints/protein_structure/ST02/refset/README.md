# ST02 formal RefSet

This directory is the ready-to-use formal input population for [ST02](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete promoted population. Do not substitute a historical small/S/4K set when a final large/L/16K population exists.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `attrition.tsv` | 384729 | columns: `entity_type`, `entity_id`, `status`, `reason` |
| `coordinates/` | 40988456 | 509 released PDB coordinate assets; every file is referenced by one formal row in `objects.tsv`. The other 1,491 source-materialization coordinates are not part of the public evaluation RefSet. |
| `dependency.tsv` | 1818 | columns: `object_a`, `object_b`, `dependency_reason`, `dependency_cluster_id` |
| `MATERIALIZATION_SUMMARY.json` | 1185 | JSON object; top-level: `artifact_id`, `edge_attrition_counts`, `formal_architecture_count`, `formal_dependency_cluster_count`, `formal_directed_endpoint_count`, `formal_object_count`, `formal_relation_count`, `formal_strict_positive_edge_count`, `maximum_negative_reuse`, `maximum_true_target_reuse`, `minimum_distinct_dependency_clusters_per_source_hsf_topology_stratum`, `necessary_same_architecture_pair_evidence_reused_count`, `parent_resolved_admitted_object_count`, `parent_unresolved_or_ambiguous_excluded_count`, `positive_candidate_edge_count`, `primary_hsf_count`, `primary_positive_topology_count`, `source_attrition_counts`, `source_hsf_topology_stratum_count`, `source_prefix_object_count`, `status` |
| `objects.tsv` | 242115 | columns: `record_id`, `source_parent_id`, `source_dependency_cluster_id`, `architecture_id`, `topology_id`, `homologous_superfamily_id`, `declared_sequence_length`, `sequence_sha256`, `coordinate_locator`, `coordinate_sha256`, `radius_of_gyration_angstrom`, `helix_fraction`, `strand_fraction`, `coil_fraction`, `parent_mapping_status` |
| `relations.tsv` | 1697489 | columns: `contrast_id`, `positive_edge_id`, `direction_id`, `anchor_object_id`, `true_target_object_id`, `negative_object_id`, `architecture_id`, `topology_id`, `source_homologous_superfamily_id`, `target_homologous_superfamily_id`, `negative_topology_id`, `source_dependency_cluster_id`, `sequence_nuisance_status`, `nuisance_match_status`, `reuse_status`, `provenance_id`, `negative_rank`, `matching_distance` |

## Handoff to model execution

Encode every structure listed in `objects.tsv`, using the PDB path in
`coordinate_locator`. Preserve the row's `record_id` as the representation ID.
The 509 released structures are exactly the dependency closure of the frozen
4,584 relations; no candidate/source-only structure needs to be encoded.

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline.
