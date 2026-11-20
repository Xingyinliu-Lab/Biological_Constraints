"""Level-0 raw-sequence ViennaRNA pairing-propensity comparator for C-RNA-05.

Input records contain only ``sequence_id`` and one raw RNA sequence. For every
sequence, ViennaRNA's unconstrained partition function yields base-pair
probabilities p_ij and this module emits q_i = sum_j p_ij. It accepts no
paired/unpaired label, coordinate map, structural constraint, chemical probing
signal, or other reference truth. Compatibility/OOS assignment and frozen
RNA05 block evaluation occur downstream of this module.
"""

import argparse
import hashlib
import json
from pathlib import Path


COMPARATOR_ID = "CMP-RNA05-VIENNA-PARTITION-PROPENSITY"
FROZEN_PARAMETERS = {
    "engine": "ViennaRNA Python bindings",
    "temperature_celsius": 37.0,
    "dangles": 2,
    "partition_function": True,
    "constraints": "FORBIDDEN",
    "chemical_reactivity": "FORBIDDEN",
    "input_alphabet": "canonical RNA A/C/G/U",
}


def _canonical_rna(sequence: str) -> str:
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("RNA sequence must be nonempty")
    if sequence != sequence.upper():
        raise ValueError("ViennaRNA comparator requires uppercase RNA input")
    if not set(sequence).issubset({"A", "C", "G", "U"}):
        raise ValueError("ViennaRNA comparator requires canonical uppercase A/C/G/U")
    return sequence


def pairing_propensity(sequence: str) -> list[float]:
    """Return q_i=sum_{j!=i} p_ij from an unconstrained partition function."""
    sequence = _canonical_rna(sequence)
    try:
        import RNA
    except ImportError as exc:
        raise RuntimeError("ViennaRNA Python bindings are required for RNA05 comparator execution") from exc
    md = RNA.md()
    md.temperature = FROZEN_PARAMETERS["temperature_celsius"]
    md.dangles = FROZEN_PARAMETERS["dangles"]
    fc = RNA.fold_compound(sequence, md)
    fc.pf()
    bpp = fc.bpp()
    n = len(sequence)
    q = []
    for i in range(1, n + 1):
        q.append(float(sum(bpp[i][j] if i < j else bpp[j][i] for j in range(1, n + 1) if j != i)))
    return q


def score_payload(payload: dict) -> dict:
    records = payload.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("input must contain a nonempty records list")
    output, seen = [], set()
    for record in records:
        sequence_id = record.get("sequence_id")
        if not isinstance(sequence_id, str) or not sequence_id or sequence_id in seen:
            raise ValueError("sequence_id must be nonempty and unique")
        seen.add(sequence_id)
        sequence = _canonical_rna(record.get("sequence"))
        output.append({
            "sequence_id": sequence_id,
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "pairing_propensity": pairing_propensity(sequence),
        })
    return {"comparator_id": COMPARATOR_ID, "parameters": FROZEN_PARAMETERS, "records": output}


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute RNA05 ViennaRNA pairing propensities")
    parser.add_argument("input", type=Path, help="JSON records: sequence_id and sequence only")
    parser.add_argument("output", type=Path, help="JSON q_i score ledger; evaluator joins truth separately")
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    args.output.write_text(json.dumps(score_payload(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
