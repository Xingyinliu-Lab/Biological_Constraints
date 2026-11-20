"""Run the public frozen DNA06 codon-local scoring protocol."""
from __future__ import annotations

import argparse
from pathlib import Path

from dna06_scoring_core import score

CONSTRAINT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFSET = CONSTRAINT_ROOT / "refset"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refset-dir", default=DEFAULT_REFSET, type=Path)
    parser.add_argument("--embeddings", required=True, type=Path)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    score(
        refset_dir=args.refset_dir,
        embeddings_path=args.embeddings,
        model_id=args.model_id,
        output_dir=args.output_dir,
        evaluation_scope="FORMAL_PROTOCOL_EVALUATION",
    )


if __name__ == "__main__":
    main()
