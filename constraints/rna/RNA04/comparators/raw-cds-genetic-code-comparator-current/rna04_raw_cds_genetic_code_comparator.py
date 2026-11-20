"""Raw-CDS Level-0 comparator for C-RNA-04.

The comparator accepts only endpoint identity plus complete WT and ALT
coding-oriented CDS strings. Its globally fixed assumptions are frame starts
at CDS nucleotide 0 and NCBI translation table 1. It finds the unique SNV,
derives its codon from that frame, and returns a score of 1 for an amino-acid
change or 0 for no amino-acid change.

It does not accept codon annotation, amino-acid annotation, mutation class, or
any experimental truth. The frozen RNA04 evaluator owns labels, complete-site
bundles, aggregation, and confidence intervals.
"""

import argparse
import json
from pathlib import Path


COMPARATOR_ID = "CMP-RNA04-RAW-CDS-GENETIC-CODE"
TRANSLATION_TABLE = 1
CANONICAL_CDS = frozenset("ACGTU")

STANDARD_CODE = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*", "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}


def _as_rna(sequence: str) -> str:
    if not isinstance(sequence, str) or not sequence or len(sequence) % 3:
        raise ValueError("CDS must be a nonempty complete codon-length string")
    sequence = sequence.upper()
    if not set(sequence).issubset(CANONICAL_CDS):
        raise ValueError("CDS contains noncanonical nucleotide symbols")
    if "T" in sequence and "U" in sequence:
        raise ValueError("CDS cannot mix T and U spellings")
    return sequence.replace("T", "U")


def _translate(codon: str) -> str:
    try:
        return STANDARD_CODE[codon]
    except KeyError as exc:
        raise ValueError(f"invalid coding codon: {codon!r}") from exc


def score_variant(wt_cds: str, alt_cds: str) -> dict:
    """Derive one endpoint score from raw WT/ALT CDS sequences only."""
    wt, alt = _as_rna(wt_cds), _as_rna(alt_cds)
    if len(wt) != len(alt):
        raise ValueError("WT and ALT CDS lengths differ")
    mismatches = [index for index, (left, right) in enumerate(zip(wt, alt)) if left != right]
    if len(mismatches) != 1:
        raise ValueError("RNA04 raw-CDS comparator requires exactly one SNV")
    position0 = mismatches[0]
    codon_start = position0 - position0 % 3
    wt_codon, alt_codon = wt[codon_start:codon_start + 3], alt[codon_start:codon_start + 3]
    wt_aa, alt_aa = _translate(wt_codon), _translate(alt_codon)
    if "*" in (wt_aa, alt_aa):
        raise ValueError("RNA04 ordinary primary excludes stop-codon consequences")
    return {
        "score": 1 if wt_aa != alt_aa else 0,
        "mismatch_position0": position0,
        "derived_wt_codon": wt_codon,
        "derived_alt_codon": alt_codon,
    }


def score_payload(payload: dict) -> dict:
    endpoints = payload.get("endpoints")
    if not isinstance(endpoints, list) or not endpoints:
        raise ValueError("input must contain a nonempty endpoints list")
    output, endpoint_ids = [], set()
    for endpoint in endpoints:
        endpoint_id = endpoint.get("endpoint_id")
        if not isinstance(endpoint_id, str) or not endpoint_id or endpoint_id in endpoint_ids:
            raise ValueError("endpoint_id must be nonempty and unique")
        endpoint_ids.add(endpoint_id)
        result = score_variant(endpoint.get("wt_cds"), endpoint.get("alt_cds"))
        output.append({"endpoint_id": endpoint_id, **result})
    return {
        "comparator_id": COMPARATOR_ID,
        "translation_table": TRANSLATION_TABLE,
        "frame_origin": "coding-oriented CDS nucleotide 0",
        "records": output,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Score raw RNA04 CDS endpoints with genetic code table 1")
    parser.add_argument("input", type=Path, help="JSON endpoint list containing only endpoint_id, wt_cds, alt_cds")
    parser.add_argument("output", type=Path, help="JSON endpoint score ledger for the frozen RNA04 evaluator")
    args = parser.parse_args()
    args.output.write_text(json.dumps(score_payload(json.loads(args.input.read_text(encoding="utf-8")),), indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
