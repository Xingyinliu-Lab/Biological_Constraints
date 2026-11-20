# DNA06 / nt-v2

This case records the model-specific Forward implementation and calls the unique public DNA06 scorer. Model environment and checkpoint: `paper_cases/model_resources/dna/nt-v2/`.

The current Forward implementation writes the formal codon-local contract directly as `row_ids + embeddings`. Canonical public-scorer outputs are stored in `output/`.
