#!/usr/bin/env bash
set -uo pipefail
root="${PS02B_WORKDIR:?set PS02B_WORKDIR}"
py="${PYTHON_BIN:-python}"
for model in esm2 prott5 esmc; do
 printf '%s\tSTARTED\t%s\n' "$model" "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
 if CUDA_VISIBLE_DEVICES=0 "$py" "$root/ps02b_forward.py" "$model" > "$root/$model.log" 2>&1; then
   printf '%s\tPASS\t%s\n' "$model" "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
 else
   printf '%s\tFAILED\t%s\n' "$model" "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
 fi
done
if "$py" "$root/ps02b_forward.py" score > "$root/score.log" 2>&1; then
 printf 'SCORE\tPASS\t%s\n' "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
else
 printf 'SCORE\tFAILED\t%s\n' "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
fi
printf 'ALL\tFINISHED\t%s\n' "$(date -u +%FT%TZ)" >> "$root/forward_status.tsv"
rm -f "$root/forward.pid"
