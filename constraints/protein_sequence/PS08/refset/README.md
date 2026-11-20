# PS08 formal RefSet

This directory is the ready-to-use formal input population for [PS08](../README.md). It is a released scientific asset, not a construction workspace or example subset.

## Use rules

- Encode the complete promoted population of 4,096 ladders and 12,288 endpoint sequences; no smaller public PS08 scoring population is authorized.
- Preserve canonical object, relation, fold and dependency identifiers exactly.
- Do not filter objects using model length, tokenizer behavior, embedding quality or model outputs.
- Read biological sequences/structures from the payload named by the object/membership tables.
- Use the owning constraint README for the scorer command and representation format.

## Files

| File | Bytes | Format and role |
|---|---:|---|
| `ENDPOINT_SEQUENCES.faa.gz` | 1144921 | FASTA; record identifier is the token following `>` and must match the relation/object tables |
| `PLANNED_POPULATION_MANIFEST.json` | 2996663 | JSON object; top-level: `canonical_members`, `kind`, `population_identity_sha256`, `status` |

## Handoff to model execution

## Handoff to scoring

After representation extraction, follow [`../scoring/README.md`](../scoring/README.md). All representations use the same public scoring pipeline. Paper-specific checkpoint and Forward settings do not alter the score calculation.
