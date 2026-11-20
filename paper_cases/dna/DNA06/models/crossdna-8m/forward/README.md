# DNA06 CrossDNA-8M Forward

```bash
python paper_cases/dna/DNA06/models/dnabert-2/forward/code/build_dna06_forward_inputs.py --refset-dir constraints/dna/DNA06/refset --output-dir work/DNA06/input
python paper_cases/dna/DNA06/models/crossdna-8m/forward/code/encode_dna06_crossdna_v011.py --model-dir <crossdna_checkpoint> --membership work/DNA06/input/DNA06_FORWARD_MEMBERSHIP.tsv --fasta work/DNA06/input/DNA06_FORWARD_INPUT.fasta --output work/DNA06/crossdna-8m/forward
python constraints/dna/DNA06/scoring/score_dna06.py --embeddings work/DNA06/crossdna-8m/forward/EMBEDDINGS.npz --model-id crossdna-8m --output-dir work/DNA06/crossdna-8m/score
```

Run these commands from the repository root. Prepare the shared input once; the builder requires `work/DNA06/input` not to exist and every model reuses it.

The adapter uses the released codon coordinates, the model's base encoding, frozen window stitching and codon-local pooling. Retain model/checkpoint identity, encoding, span mapping, windowing and pooling as Forward evidence.
