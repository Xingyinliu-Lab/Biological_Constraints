# DNA06 NT-v2 Forward

```bash
python paper_cases/dna/DNA06/models/dnabert-2/forward/code/build_dna06_forward_inputs.py --refset-dir constraints/dna/DNA06/refset --output-dir work/DNA06/input
DNA06_MODEL_PATH=<nt_v2_checkpoint> python paper_cases/dna/DNA06/models/nt-v2/forward/code/encode_dna06_transformer_models_v011.py --model nt --membership work/DNA06/input/DNA06_FORWARD_MEMBERSHIP.tsv --fasta work/DNA06/input/DNA06_FORWARD_INPUT.fasta --output work/DNA06/nt-v2/forward
python constraints/dna/DNA06/scoring/score_dna06.py --embeddings work/DNA06/nt-v2/forward/EMBEDDINGS.npz --model-id nt-v2 --output-dir work/DNA06/nt-v2/score
```

Run these commands from the repository root. Prepare the shared input once; the builder requires `work/DNA06/input` not to exist and every model reuses it.

The adapter uses the released codon coordinates, model token spans, frozen window stitching and codon-local pooling. Retain model/checkpoint/tokenizer identity, span mapping, windowing and pooling as Forward evidence.
