#!/usr/bin/env python3
"""Run the frozen RNA06 ViennaRNA partner-probability comparator."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import random
from pathlib import Path


HERE = Path(__file__).resolve().parent
CONSTRAINT = HERE.parents[1]
PUBLIC_SCORER = CONSTRAINT / "scoring" / "score_rna06.py"


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha_object(value):
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def load_public_scorer():
    spec = importlib.util.spec_from_file_location("rna06_public_scorer", PUBLIC_SCORER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def aggregate_with_public_scorer(records, public_scorer):
    """Adapt ViennaRNA measurements to the one public RNA06 scoring hierarchy."""
    rows = []
    scores = []
    true_score_by_endpoint = {}
    for record in records:
        member = record["member"]
        endpoint = member["directed_endpoint_id"]
        true_score = float(record["p_anchor_true_partner"])
        is_new_endpoint = endpoint not in true_score_by_endpoint
        previous = true_score_by_endpoint.setdefault(endpoint, true_score)
        if previous != true_score:
            raise ValueError(f"true-partner score differs within endpoint: {endpoint}")
        common = {
            "outer_dependency_cluster_id": member["dependency_cluster_id"],
            "parent_id": member["parent_id"],
            "edge_bundle_id": member["true_edge_id"],
            "directed_endpoint_id": endpoint,
        }
        if is_new_endpoint:
            rows.append({**common, "class": "true_partner"})
            scores.append(true_score)
        rows.append({**common, "class": "paired_nonpartner"})
        scores.append(float(record["p_anchor_nonpartner"]))
    return public_scorer.aggregate(rows, scores)


def read_fasta_by_sequence_hash(path: Path, sequence_sha256):
    records = {}
    header = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                sequence = "".join(chunks).upper()
                records[sequence_sha256(sequence)] = sequence
            header = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        sequence = "".join(chunks).upper()
        records[sequence_sha256(sequence)] = sequence
    return records


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the frozen RNA06 ViennaRNA comparator from released inputs."
    )
    parser.add_argument("--sequences", type=Path, default=CONSTRAINT / "refset" / "rna_sequences.fasta")
    parser.add_argument("--objects", type=Path, default=CONSTRAINT / "refset" / "rna_objects.tsv")
    parser.add_argument(
        "--pair-state-maps", type=Path, default=CONSTRAINT / "refset" / "RNA06_FROZEN_PAIR_STATE_MAPS.tsv"
    )
    parser.add_argument(
        "--directed-endpoints",
        type=Path,
        default=CONSTRAINT / "refset" / "RNA06_FORMAL_DIRECTED_ENDPOINT_LEDGER.tsv",
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=CONSTRAINT / "refset" / "RNA06_FORMAL_PAIRED_NONPARTNER_CANDIDATE_LEDGER.tsv",
    )
    parser.add_argument(
        "--edges",
        type=Path,
        default=CONSTRAINT / "refset" / "RNA06_FORMAL_RECIPROCAL_EDGE_BUNDLE_LEDGER.tsv",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    observable_spec = importlib.util.spec_from_file_location(
        "rna06_vienna_observable", HERE / "vienna_partition_partner_probability.py"
    )
    observable = importlib.util.module_from_spec(observable_spec)
    observable_spec.loader.exec_module(observable)
    public_scorer = load_public_scorer()

    objects = {row["record_id"]: row["sequence_sha256"] for row in read_tsv(args.objects)}
    sequences = read_fasta_by_sequence_hash(args.sequences, observable.sequence_sha256)
    pairsets = {}
    for row in read_tsv(args.pair_state_maps):
        pairs = [(int(x[0]), int(x[1]), str(x[2])) for x in json.loads(row["v1_eligible_positive_pairs_json"])]
        pairsets[(row["record_id"], row["structure_map_id"], row["condition_id"])] = pairs

    directed = {row["directed_endpoint_id"]: row for row in read_tsv(args.directed_endpoints)}
    candidates = read_tsv(args.candidates)
    edges = {row["edge_bundle_id"]: row for row in read_tsv(args.edges)}
    edge_oos = []
    compatible_edges = set()
    for edge_id, edge in edges.items():
        key = (edge["record_id"], edge["primary_structure_map_id"], edge["primary_condition_id"])
        truth = pairsets.get(key)
        sequence_hash = objects.get(edge["record_id"])
        if truth is None or sequence_hash not in sequences:
            raise ValueError(f"missing frozen source object {key!r}")
        a, b = sorted((int(edge["position_i0"]), int(edge["position_j0"])))
        match = [x for x in truth if frozenset((x[0], x[1])) == frozenset((a, b))]
        reasons = []
        if len(match) != 1:
            reasons.append("TRUE_EDGE_NOT_IN_V1_PAIR_UNIVERSE")
        else:
            pair_type = match[0][2]
            chemistry = sequences[sequence_hash][a] + sequences[sequence_hash][b]
            if pair_type != "cWW":
                reasons.append("UNSUPPORTED_NONCANONICAL_OR_TERTIARY_PAIR_TYPE")
            if chemistry not in {"AU", "UA", "GC", "CG", "GU", "UG"}:
                reasons.append("UNSUPPORTED_PAIR_CHEMISTRY")
            for x, y, _ in truth:
                c, d = sorted((x, y))
                if len({a, b, c, d}) == 4 and ((a < c < b < d) or (c < a < d < b)):
                    reasons.append("PSEUDOKNOT_CROSSING_TRUE_EDGE")
                    break
        if reasons:
            edge_oos.append({"edge_bundle_id": edge_id, "reasons": reasons})
        else:
            compatible_edges.add(edge_id)

    cache = {}
    sequence_audit = []
    for record_id in sorted({row["record_id"] for row in directed.values()}):
        sequence_hash = objects.get(record_id)
        if sequence_hash is None or sequence_hash not in sequences:
            raise ValueError(f"missing raw sequence {record_id}")
        bpp, length = observable.partition_probabilities(sequences[sequence_hash])
        cache[record_id] = (bpp, length)
        sequence_audit.append(
            {"record_id": record_id, "sequence_sha256": sequence_hash, "length_nt": length}
        )

    import RNA

    RNA.params_load_RNA_Turner2004()
    parameter_artifact = args.output / "RNA_Turner2004_ViennaRNA_2.4.17.par"
    if RNA.params_save(str(parameter_artifact)) != 1:
        raise RuntimeError("failed to materialize ViennaRNA parameter artifact")
    parameter_sha256 = hashlib.sha256(parameter_artifact.read_bytes()).hexdigest()

    records = []
    for candidate in candidates:
        if candidate["edge_bundle_id"] not in compatible_edges:
            continue
        endpoint = directed[candidate["directed_endpoint_id"]]
        bpp, length = cache[candidate["record_id"]]
        p_true = observable.partner_probability(
            bpp, length, int(candidate["anchor_position0"]), int(candidate["true_partner_position0"])
        )
        p_nonpartner = observable.partner_probability(
            bpp, length, int(candidate["anchor_position0"]), int(candidate["paired_nonpartner_position0"])
        )
        credit = 1.0 if p_true > p_nonpartner else (-1.0 if p_true < p_nonpartner else 0.0)
        member = {
            "leaf_id": candidate["candidate_id"],
            "matched_nonpartner_id": candidate["candidate_id"],
            "directed_endpoint_id": candidate["directed_endpoint_id"],
            "true_edge_id": candidate["edge_bundle_id"],
            "parent_id": endpoint["parent_id"],
            "dependency_cluster_id": endpoint["outer_dependency_cluster_id"],
        }
        records.append(
            {
                "leaf_id": candidate["candidate_id"],
                "member": member,
                "F_credit": credit,
                "p_anchor_true_partner": p_true,
                "p_anchor_nonpartner": p_nonpartner,
            }
        )
    if not records:
        raise ValueError("no comparator-compatible candidates")

    population_identity = sha_object(sorted([record["member"] for record in records], key=canonical_json))
    (
        score,
        cluster_scores,
        endpoint_scores,
        endpoint_raw,
        edge_scores,
        parent_scores,
        parent_cluster,
    ) = aggregate_with_public_scorer(
        records, public_scorer
    )
    # Bootstrap ordering is part of the frozen deterministic procedure: the
    # seeded draws are indexed by canonical dependency-cluster ID order.
    cluster_values = [cluster_scores[key] for key in sorted(cluster_scores)]
    rng = random.Random(20261004)
    draws = sorted(
        math.fsum(rng.choice(cluster_values) for _ in cluster_values) / len(cluster_values)
        for _ in range(20000)
    )
    ci = [draws[int(0.025 * (len(draws) - 1))], draws[int(0.975 * (len(draws) - 1))]]

    audit = {
        "comparator_id": "CMP-RNA06-VIENNA-PARTITION-PARTNER-PROBABILITY",
        "input_boundary": "raw RNA sequence only; truth and constraints are not supplied to ViennaRNA",
        "formal_edges_total": len(edges),
        "compatible_edges": len(compatible_edges),
        "oos_edges": len(edge_oos),
        "oos": edge_oos,
        "sequence_audit": sequence_audit,
        "thermodynamic_parameter_set": "RNA_Turner2004",
        "thermodynamic_parameter_artifact_sha256": parameter_sha256,
    }
    result = {
        "comparator_id": audit["comparator_id"],
        "engine": {
            "name": "ViennaRNA Python bindings",
            "version": RNA.__version__,
            "temperature_c": 37.0,
            "dangles": 2,
            "thermodynamic_parameter_set": "RNA_Turner2004",
            "thermodynamic_parameter_artifact_sha256": parameter_sha256,
            "partition_function": True,
            "constraints": "FORBIDDEN",
            "chemical_reactivity": "FORBIDDEN",
        },
        "population": {
            "formal_edge_count": len(edges),
            "compatible_edge_count": len(compatible_edges),
            "oos_edge_count": len(edge_oos),
            "scored_candidate_contrasts": len(records),
            "scored_directed_endpoint_count": len({r["member"]["directed_endpoint_id"] for r in records}),
            "scored_parent_count": len({r["member"]["parent_id"] for r in records}),
            "scored_outer_cluster_count": len(cluster_values),
            "compatible_population_identity_sha256": population_identity,
        },
        "result": {
            "F": score,
            "F_95_CI": ci,
            "bootstrap_unit": "outer_dependency_cluster_F",
            "bootstrap_replicates": 20000,
            "bootstrap_seed": 20261004,
        },
        "execution_identity": {
            "comparator_code_sha256": hashlib.sha256((HERE / "vienna_partition_partner_probability.py").read_bytes()).hexdigest(),
            "public_scorer": str(PUBLIC_SCORER.relative_to(CONSTRAINT.parents[2])).replace("\\", "/"),
            "public_scorer_sha256": hashlib.sha256(PUBLIC_SCORER.read_bytes()).hexdigest(),
            "aggregation_output_sha256": sha_object(
                {
                    "cluster_scores": cluster_scores,
                    "endpoint_scores": endpoint_scores,
                    "endpoint_raw": endpoint_raw,
                    "edge_scores": edge_scores,
                    "parent_scores": parent_scores,
                    "parent_cluster": parent_cluster,
                }
            ),
        },
    }
    (args.output / "RNA06_VIENNA_COMPATIBILITY_AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n")
    (args.output / "RNA06_VIENNA_SCORE_LEDGER.json").write_text(
        json.dumps({"records": records}, separators=(",", ":")) + "\n"
    )
    candidate_rows=[]
    for record in records:
        member=record['member']; true=float(record['p_anchor_true_partner']); nonpartner=float(record['p_anchor_nonpartner'])
        candidate_rows.append({'candidate_id':record['leaf_id'],'directed_endpoint_id':member['directed_endpoint_id'],'edge_bundle_id':member['true_edge_id'],'parent_id':member['parent_id'],'outer_dependency_cluster_id':member['dependency_cluster_id'],'true_score':true,'nonpartner_score':nonpartner,'margin':true-nonpartner,'F_credit':float(record['F_credit'])})
    public_scorer.write_tsv(args.output/'CANDIDATE_SCORES.tsv',candidate_rows,list(candidate_rows[0]))
    endpoint_meta={record['member']['directed_endpoint_id']:record['member'] for record in records}
    endpoint_rows=[{'directed_endpoint_id':endpoint,'edge_bundle_id':endpoint_meta[endpoint]['true_edge_id'],'parent_id':endpoint_meta[endpoint]['parent_id'],'outer_dependency_cluster_id':endpoint_meta[endpoint]['dependency_cluster_id'],**endpoint_raw[endpoint],'individual_F':value} for endpoint,value in sorted(endpoint_scores.items())]
    public_scorer.write_tsv(args.output/'INDIVIDUAL_SCORES.tsv',endpoint_rows,list(endpoint_rows[0]))
    edge_meta={record['member']['true_edge_id']:record['member'] for record in records}
    edge_rows=[{'edge_bundle_id':edge,'parent_id':edge_meta[edge]['parent_id'],'outer_dependency_cluster_id':edge_meta[edge]['dependency_cluster_id'],'edge_F':value} for edge,value in sorted(edge_scores.items())]
    public_scorer.write_tsv(args.output/'EDGE_SCORE_LEDGER.tsv',edge_rows,list(edge_rows[0]))
    parent_rows=[{'parent_id':parent,'outer_dependency_cluster_id':parent_cluster[parent],'parent_F':value} for parent,value in sorted(parent_scores.items())]
    public_scorer.write_tsv(args.output/'PARENT_SCORE_LEDGER.tsv',parent_rows,list(parent_rows[0]))
    cluster_rows=[{'outer_dependency_cluster_id':cluster,'cluster_F':value} for cluster,value in sorted(cluster_scores.items())]
    public_scorer.write_tsv(args.output/'CLUSTER_SCORE_LEDGER.tsv',cluster_rows,list(cluster_rows[0]))
    result["formal_interface"] = "RNA06_F_ONLY_v0.1.3"
    result["score_transform_version"] = "direct_F_credit_v1"
    (args.output / "FORMAL_RESULT.json").write_text(json.dumps(result, indent=2) + "\n")
    files=[path for path in args.output.iterdir() if path.is_file() and path.name!='OUTPUT_MANIFEST.json']
    (args.output/'OUTPUT_MANIFEST.json').write_text(json.dumps({'formal_interface':'RNA06_F_ONLY_v0.1.3','files':[{'path':path.name,'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()} for path in sorted(files)]},indent=2)+'\n')
    print(json.dumps({"F": score, "F_95_CI": ci, "records": len(records)}, sort_keys=True))


if __name__ == "__main__":
    main()
