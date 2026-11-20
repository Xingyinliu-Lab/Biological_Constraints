# Forward

The formal PS10 sequence is `normalized_peptide_identity`. Prepare the ESM3 FASTA, execute the retained ESM3 extractor, then export its Safetensors/index pair to the public NPZ contract:

```bash
python paper_cases/protein_sequence/PS10/models/esm3-seq/forward/code/ps10_esm3_handoff.py prepare-fasta \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --output work/PS10/esm3-seq/forward/PS10_PEPTIDES.faa

python paper_cases/protein_sequence/PS10/models/esm3-seq/forward/code/esm3_seq_mean_extract.py \
  --fasta work/PS10/esm3-seq/forward/PS10_PEPTIDES.faa \
  --id-column record_id \
  --output-dir work/PS10/esm3-seq/forward/raw

python paper_cases/protein_sequence/PS10/models/esm3-seq/forward/code/ps10_esm3_handoff.py export-npz \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --embeddings work/PS10/esm3-seq/forward/raw/embeddings.safetensors \
  --index work/PS10/esm3-seq/forward/raw/embedding_index.tsv \
  --output work/PS10/esm3-seq/forward/embeddings.npz
```

The NPZ contains exactly `record_id` and `embedding`, preserves all 2,031 record identities, and is the input used by the public scoring command.
