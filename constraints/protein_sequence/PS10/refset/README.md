# PS10 formal RefSet

This directory contains the complete ready-to-Forward PS10 population and its frozen study-aware assignments. It contains no historical candidate pool.

| File | Role |
|---|---|
| `PEPTIDE_REFSET_RECORDS.tsv` | Unique peptide records, canonical `normalized_peptide_identity` sequences, experimental classes and source-study units. These are the model Forward inputs. |
| `PRIMARY_ASSIGNMENTS.tsv` | Frozen train/test/purged role for every record in each of 22 study-held-out folds; defines the registered result. |
| `SENSITIVITY_ASSIGNMENTS.tsv` | Separate frozen sensitivity design; never substitutes for the primary result. |

Encode the `normalized_peptide_identity` value for every unique record and return one sequence-level vector per `record_id`. Its `length` and `peptide_identity_sha256` columns identify the exact formal sequence. The public scorer consumes the assignments itself; users must not manufacture fold-specific predictions or choose scoring parameters.
