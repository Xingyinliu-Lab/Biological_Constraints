# ST05 paper implementation cases

Five registered structure-model Forward implementations submit the same four-array representation contract to the unique public ST05 invariant scorer. Comparator execution and results remain only in the Constraint release. Large representation NPZ files are not duplicated in Git.

| Model | F | F 95% CI | Public parity |
|---|---:|---|---|
| esm-if1 | 0.4364406779661017 | [0.4311440677966102, 0.4406779661016949] | PASS |
| esm3-struct | 0.2801906779661017 | [0.2722457627118644, 0.288135593220339] | PASS |
| proteinmpnn-ca | 0.5540254237288136 | [0.5074152542372882, 0.6016949152542372] | PASS |
| protoken | 0.1186440677966102 | [0.0858050847457627, 0.16101694915254236] | PASS |
| sslm | 0.954978813559322 | [0.9322033898305084, 0.9751059322033898] | PASS |

Each case contains the actual Forward scripts, the exact scoring command, and published F-only outputs. Historical embedding NPZ files and their source hashes are not distributed. Public reproduction begins with the released RefSet and model Forward manual, then runs the public scorer to regenerate input identities, the ID-level F ledger and the aggregate result; the published summaries do not claim byte-level identity to non-public historical NPZ files.
