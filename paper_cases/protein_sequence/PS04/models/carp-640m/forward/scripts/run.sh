#!/usr/bin/env bash
set -euo pipefail
root="${PS04_WORKDIR:?set PS04_WORKDIR}"
py="${PYTHON_BIN:-python}"
export TORCH_HOME="${TORCH_HOME:?set TORCH_HOME}"
printf "carp\tSTARTED\t%s\n" "$(date -u +%FT%TZ)" > "$root/status.tsv"
CUDA_VISIBLE_DEVICES=0 "$py" "$root/ps04_carp640_forward.py" carp > "$root/carp.log" 2>&1
"$py" "$root/ps04_carp640_forward.py" score > "$root/score.log" 2>&1
printf "ALL_FINISHED\t%s\n" "$(date -u +%FT%TZ)" >> "$root/status.tsv"
rm -f "$root/forward.pid"
