# DNA02-02 / crossdna-8m execution commands

Run these commands from the repository root in the environment documented by `paper_cases/model_resources/dna/crossdna-8m/`.

First create the exact three-column table expected by the archived Forward implementation:

```bash
python paper_cases/dna/DNA02-02/models/crossdna-8m/forward/scripts/build_forward_object_table.py \
  --sequence-objects constraints/dna/DNA02-02/refset/SEQUENCE_OBJECT.tsv \
  --output work/DNA02-02/crossdna-8m/FORWARD_OBJECT_TABLE.tsv
```

The command below reproduces the historical execution policy: it loads the checkpoint once and performs two eval-mode Forward passes without reloading between them. Adding `--reload-between-passes` defines a separate optional robustness run and does not reproduce the archived paper representations:

```bash
python paper_cases/dna/DNA02-02/models/crossdna-8m/forward/scripts/encode_dna01_crossdna.py \
  --object-table work/DNA02-02/crossdna-8m/FORWARD_OBJECT_TABLE.tsv \
  --fasta constraints/dna/DNA02-02/refset/DNA02_02_SEQUENCE_OBJECTS.fasta \
  --model-dir <CROSSDNA_MODEL_DIR> \
  --output work/DNA02-02/crossdna-8m/forward
```

Use the model-specific script in `forward/scripts/` and the checkpoint paths described by the model resource. The Forward output directory must contain `EMBEDDINGS.npz` and `REPEAT_EMBEDDINGS.npz`.

Then run the unique public scorer:

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings work/DNA02-02/crossdna-8m/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA02-02/crossdna-8m/forward/REPEAT_EMBEDDINGS.npz \
  --model-id crossdna-8m \
  --output work/DNA02-02/crossdna-8m/score
```
