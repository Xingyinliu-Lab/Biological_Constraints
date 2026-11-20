# SSLM Forward

First derive one Foldseek 3Di record per released coordinate.

```bash
mkdir -p work/ST01/sslm/3di
foldseek createdb constraints/protein_structure/ST01/refset/coordinates work/ST01/sslm/3di/db
foldseek lndb work/ST01/sslm/3di/db_h work/ST01/sslm/3di/db_ss_h
foldseek convert2fasta work/ST01/sslm/3di/db_ss work/ST01/sslm/3di/structures_3di.fasta
python paper_cases/protein_structure/ST01/models/sslm/forward/scripts/sslm_mean_extract.py --fasta work/ST01/sslm/3di/structures_3di.fasta --objects constraints/protein_structure/ST01/refset/objects.tsv --model-dir /path/to/sslm --output work/ST01/sslm/embeddings.npz
```

The extractor maps the Foldseek filename stem to the canonical ST01 ID and mean-pools valid 3Di token states. Score with [`../scoring_call/COMMANDS.md`](../scoring_call/COMMANDS.md). The model environment is documented in [`model_resources/protein_structure/sslm`](../../../../../model_resources/protein_structure/sslm/README.md).
