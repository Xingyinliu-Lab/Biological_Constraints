#!/usr/bin/env bash
set -euo pipefail
exec python paper_cases/protein_sequence/PS01/shared_forward/ps01_forward_public.py --model prott5 "$@"
