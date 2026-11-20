"""Exact-Hamming Level-0 oracle for C-RNA-01-02.

The input JSON contains a nonempty ``leaves`` list. Every leaf must contain
the frozen ladder provenance fields, raw canonical RNA sequences ``x0/x1/x2``,
and SHA256 identities ``x0_hash/x1_hash/x2_hash``. The emitted distances are
integer Hamming counts and remain compatible with the generic evaluator's
``d01/d12/d02`` fields.

Each input must be an exact single-edit ladder, H01/H12/H02 = 1/1/2. Score
semantics are permanently bound to ``oracle_tau_delta = 0``; callers must not
substitute a model-calibrated tolerance.
"""

import argparse
import hashlib
import json
from pathlib import Path


COMPARATOR_ID = "CMP-RNA01-02-EXACT-HAMMING"
ORACLE_TAU_DELTA = 0
CANONICAL_RNA = frozenset("ACGU")

IDENTITY_FIELDS = (
    "leaf_id",
    "ladder_id",
    "first_edit_id",
    "second_edit_id",
    "source_parent_id",
    "outer_cluster_id",
    "tau_stratum_id",
    "x0_hash",
    "x1_hash",
    "x2_hash",
)


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def exact_hamming(left: str, right: str) -> int:
    """Return exact integer Hamming distance after canonical-RNA validation."""
    if not isinstance(left, str) or not isinstance(right, str):
        raise ValueError("RNA sequences must be strings")
    if not left or len(left) != len(right):
        raise ValueError("RNA01-02 requires equal, nonempty sequences")
    if not set(left).issubset(CANONICAL_RNA) or not set(right).issubset(CANONICAL_RNA):
        raise ValueError("RNA01-02 requires canonical uppercase A/C/G/U sequences")
    return sum(a != b for a, b in zip(left, right))


def _validated_sequence(leaf: dict, state: str) -> str:
    sequence = leaf[state]
    expected = leaf[f"{state}_hash"]
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError(f"invalid {state}_hash")
    if sequence_sha256(sequence) != expected:
        raise ValueError(f"{state}_hash does not bind the supplied {state} sequence")
    return sequence


def score_leaf(leaf: dict) -> dict:
    missing = [key for key in (*IDENTITY_FIELDS, "x0", "x1", "x2") if key not in leaf]
    if missing:
        raise ValueError("missing required field(s): " + ", ".join(missing))
    x0 = _validated_sequence(leaf, "x0")
    x1 = _validated_sequence(leaf, "x1")
    x2 = _validated_sequence(leaf, "x2")
    h01, h12, h02 = exact_hamming(x0, x1), exact_hamming(x1, x2), exact_hamming(x0, x2)
    if (h01, h12, h02) != (1, 1, 2):
        raise ValueError("invalid RNA01-02 frozen ladder: expected exact Hamming 1/1/2")
    out = {key: leaf[key] for key in IDENTITY_FIELDS}
    # d* names preserve compatibility with the frozen generic evaluator.
    out.update(h01=h01, h12=h12, h02=h02, d01=h01, d12=h12, d02=h02)
    return out


def score_payload(payload: dict) -> dict:
    leaves = payload.get("leaves")
    if not isinstance(leaves, list) or not leaves:
        raise ValueError("input must contain a nonempty leaves list")
    scored = [score_leaf(leaf) for leaf in leaves]
    for key in ("leaf_id", "ladder_id"):
        values = [row[key] for row in scored]
        if len(set(values)) != len(values):
            raise ValueError(f"duplicate {key}")
    return {
        "comparator_id": COMPARATOR_ID,
        "oracle_tau_delta": ORACLE_TAU_DELTA,
        "leaves": scored,
    }


def score_frozen_ledgers(ladder_rows: list[dict], object_rows: list[dict]) -> dict:
    """Adapt the frozen RNA01-02 ladder/object ledgers without exposing sequence output."""
    objects = {row["object_id"]: row for row in object_rows}
    if len(objects) != len(object_rows):
        raise ValueError("duplicate object_id in frozen object ledger")
    leaves = []
    for ladder in ladder_rows:
        object_ids = {state: ladder[f"{state}_object_id"] for state in ("x0", "x1", "x2")}
        if any(object_id not in objects for object_id in object_ids.values()):
            raise ValueError("ladder references an absent frozen object")
        leaf = {
            "leaf_id": ladder["leaf_id"],
            "ladder_id": ladder["leaf_id"],
            "first_edit_id": ladder["first_edit_id"],
            "second_edit_id": ladder["nested_extension_id"],
            "source_parent_id": ladder["source_parent_id"],
            "outer_cluster_id": ladder["outer_cluster_id"],
            "tau_stratum_id": ladder["tau_stratum_id"],
        }
        for state, object_id in object_ids.items():
            obj = objects[object_id]
            ledger_hash = ladder[f"{state}_sequence_sha256"]
            if obj["sequence_sha256"] != ledger_hash:
                raise ValueError(f"{state} object/ladder sequence hash disagreement")
            leaf[state] = obj["sequence"]
            leaf[f"{state}_hash"] = ledger_hash
        scored = score_leaf(leaf)
        scored.update({f"{state}_object_id": object_id for state, object_id in object_ids.items()})
        leaves.append(scored)
    for key in ("leaf_id", "ladder_id"):
        values = [row[key] for row in leaves]
        if len(set(values)) != len(values):
            raise ValueError(f"duplicate {key}")
    return {
        "comparator_id": COMPARATOR_ID,
        "oracle_tau_delta": ORACLE_TAU_DELTA,
        "leaves": leaves,
    }


def require_oracle_tau(tau_delta: int | float) -> None:
    """Reject scorer invocations not bound to the oracle tolerance."""
    if tau_delta != ORACLE_TAU_DELTA:
        raise ValueError("exact-Hamming oracle requires comparator tau_delta=0")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create RNA01-02 exact-Hamming oracle payload")
    parser.add_argument("input", type=Path, help="JSON ladder payload with x0/x1/x2 and sequence hashes")
    parser.add_argument("output", type=Path, help="JSON integer-distance payload for the RNA01-02 evaluator")
    args = parser.parse_args()
    result = score_payload(json.loads(args.input.read_text(encoding="utf-8")))
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
