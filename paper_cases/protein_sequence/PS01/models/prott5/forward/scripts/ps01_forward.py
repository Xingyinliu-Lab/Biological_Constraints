#!/usr/bin/env python3
"""Portable PS01 Forward wrapper for prott5."""
from pathlib import Path
import runpy
import sys

sys.argv[1:1] = ["--model", "prott5"]
adapter = Path(__file__).resolve().parents[4] / "shared_forward" / "ps01_forward_public.py"
runpy.run_path(str(adapter), run_name="__main__")
