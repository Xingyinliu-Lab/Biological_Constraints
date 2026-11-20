"""Raw-sequence ViennaRNA partner-probability comparator for C-RNA-06.

This constraint-local Level-0 method accepts only a raw canonical RNA sequence.
It computes an unconstrained partition function and exposes p_ij.  Frozen RNA06
truth, true-partner identities, and nonpartner candidates are joined only by the
downstream evaluator.
"""

import hashlib


COMPARATOR_ID = "CMP-RNA06-VIENNA-PARTITION-PARTNER-PROBABILITY"
FROZEN_PARAMETERS = {
    "engine": "ViennaRNA Python bindings",
    "temperature_celsius": 37.0,
    "dangles": 2,
    "thermodynamic_parameter_set": "RNA_Turner2004 built-in parameter set",
    "partition_function": True,
    "constraints": "FORBIDDEN",
    "chemical_reactivity": "FORBIDDEN",
    "input_alphabet": "canonical RNA A/C/G/U",
}


def _canonical_rna(sequence: str) -> str:
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("RNA sequence must be nonempty")
    if sequence != sequence.upper() or not set(sequence).issubset({"A", "C", "G", "U"}):
        raise ValueError("ViennaRNA comparator requires canonical uppercase A/C/G/U")
    return sequence


def partition_probabilities(sequence: str):
    """Return the ViennaRNA BPP matrix for one raw RNA sequence."""
    sequence = _canonical_rna(sequence)
    try:
        import RNA
    except ImportError as exc:
        raise RuntimeError("ViennaRNA Python bindings are required for RNA06 comparator execution") from exc
    RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    md.temperature = FROZEN_PARAMETERS["temperature_celsius"]
    md.dangles = FROZEN_PARAMETERS["dangles"]
    fc = RNA.fold_compound(sequence, md)
    fc.pf()
    return fc.bpp(), len(sequence)


def partner_probability(bpp, length: int, position_i0: int, position_j0: int) -> float:
    """Read p_ij using external zero-based RNA coordinates."""
    if not (0 <= position_i0 < length and 0 <= position_j0 < length) or position_i0 == position_j0:
        raise ValueError("invalid RNA06 partner coordinate")
    i, j = sorted((position_i0 + 1, position_j0 + 1))
    return float(bpp[i][j])


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(_canonical_rna(sequence).encode()).hexdigest()
