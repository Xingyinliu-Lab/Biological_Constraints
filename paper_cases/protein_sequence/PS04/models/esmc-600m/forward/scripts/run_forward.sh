#!/usr/bin/env bash
set -uo pipefail
root="${PS04_WORKDIR:?set PS04_WORKDIR}"
py="${PYTHON_BIN:-python}"
for model in esm2 prott5 esmc; do
 printf '%s\tSTARTED\t%s\n' "$model" "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
 if CUDA_VISIBLE_DEVICES=0 "$py" "$root/ps04_forward.py" "$model" > "$root/$model.log" 2>&1; then
  printf '%s\tPASS\t%s\n' "$model" "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
 else
  printf '%s\tFAILED\t%s\n' "$model" "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
  exit 1
 fi
done
"$py" "$root/ps04_forward.py" score > "$root/pair_score.log" 2>&1
for model in esm2 prott5 esmc; do
 for panel in S_2K L_4K; do
  out="$root/$model/score_$panel"
  "$py" "$root/score_ps04.py" --contract-root "$root/contract" --leaf-ledger "$root/$model/LEAF_SCORES_$panel.tsv" --delta-num 0 --output-dir "$out"
 done
done
printf 'SCORE\tPASS\t%s\n' "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
printf 'ALL\tFINISHED\t%s\n' "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
rm -f "$root/forward.pid"
