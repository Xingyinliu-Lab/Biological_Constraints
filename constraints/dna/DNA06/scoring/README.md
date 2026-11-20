# DNA06 scoring

`score_dna06.py` is the only public model-scoring entry. It accepts the codon-local NPZ specified by [`../README.md`](../README.md) and performs every Constraint-defined operation after model Forward.

```bash
python constraints/dna/DNA06/scoring/score_dna06.py \
  --embeddings <EMBEDDINGS_NPZ> \
  --model-id <MODEL_ID> \
  --output-dir <NEW_OUTPUT_DIRECTORY>
```

The released RefSet is located automatically. `dna06_scoring_core.py` performs input validation, feature construction, frozen cross-fitting and formal output generation; `dna06_numeric_core.py` contains the fixed numerical primitives. There is no model-specific scoring path.

The released biological targets are in [`../refset/FORWARD_TARGETS.tsv`](../refset/FORWARD_TARGETS.tsv); the full NPZ and mapping-evidence requirement is in [`INPUT_CONTRACT.json`](INPUT_CONTRACT.json). Tokenization and token-to-nucleotide mapping remain the submitter's model-specific responsibility.

The frozen ridge objective is the **sum** of per-variant binary log-loss plus
`0.5 * lambda * ||beta||^2`, with the intercept unpenalized. An earlier
`mean_binary_log_loss` label was corrected as method metadata on 2026-10-19;
the historical executed code, current public core and registered FM results all
use the summed objective.
