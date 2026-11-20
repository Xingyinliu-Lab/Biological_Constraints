# DNA02-02 / dnabert-2 execution commands

Run these commands from the repository root in the environment documented by `paper_cases/model_resources/dna/dnabert-2/`.

First create the exact three-column table expected by the archived Forward implementation:

```bash
python paper_cases/dna/DNA02-02/models/dnabert-2/forward/scripts/build_forward_object_table.py \
  --sequence-objects constraints/dna/DNA02-02/refset/SEQUENCE_OBJECT.tsv \
  --output work/DNA02-02/dnabert-2/FORWARD_OBJECT_TABLE.tsv
```

The command below reproduces the historical execution policy: it loads the checkpoint once and performs two eval-mode Forward passes without reloading between them. Adding `--reload-between-passes` defines a separate optional robustness run and does not reproduce the archived paper representations:

```bash
DNA02_02_DNABERT2_DIRECT_HF=1 python paper_cases/dna/DNA02-02/models/dnabert-2/forward/scripts/encode_dna01_transformer_models.py \
  --model dnabert2 \
  --object-table work/DNA02-02/dnabert-2/FORWARD_OBJECT_TABLE.tsv \
  --fasta constraints/dna/DNA02-02/refset/DNA02_02_SEQUENCE_OBJECTS.fasta \
  --output work/DNA02-02/dnabert-2/forward \
  --strict-repeat
```

Use the model-specific script in `forward/scripts/` and the checkpoint paths described by the model resource. The Forward output directory must contain `EMBEDDINGS.npz` and `REPEAT_EMBEDDINGS.npz`.

Then run the unique public scorer:

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings work/DNA02-02/dnabert-2/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA02-02/dnabert-2/forward/REPEAT_EMBEDDINGS.npz \
  --model-id dnabert-2 \
  --output work/DNA02-02/dnabert-2/score
```
