from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

BOOTSTRAP_DRAWS = 20_000
CALIBRATION_POLICY = "maximum absolute difference between repeat complete-ladder contrasts"
REGISTERED_POPULATION_IDENTITY_SHA256 = "686a035451271f480783291f316862bd4932506735c4d24d8f60aef69ab8b558"
CANONICAL_RNA = frozenset("ACGU")
POPULATION_MEMBER_KEYS = (
    "first_edit_id", "leaf_id", "nested_extension_id", "outer_cluster_id",
    "source_parent_id", "source_unit_id", "x0_object_id", "x0_sequence_sha256",
    "x1_object_id", "x1_sequence_sha256", "x2_object_id", "x2_sequence_sha256",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_vectors(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {"sequence_sha256", "vectors"}:
            raise ValueError(f"{path}: expected exactly sequence_sha256 and vectors")
        raw_ids = np.asarray(data["sequence_sha256"])
        raw_vectors = np.asarray(data["vectors"])
        if raw_ids.ndim != 1 or raw_ids.dtype.kind not in "SU":
            raise ValueError(f"{path}: sequence_sha256 must be a one-dimensional string array")
        if raw_vectors.dtype.kind not in "iuf":
            raise ValueError(f"{path}: vectors must use a real numeric dtype")
        ids = [str(value) for value in raw_ids]
        vectors = np.asarray(raw_vectors, dtype=np.float64)
    if vectors.ndim != 2 or vectors.shape[1] <= 0 or vectors.shape[0] != len(ids):
        raise ValueError(f"{path}: invalid vector matrix")
    if len(set(ids)) != len(ids):
        raise ValueError(f"{path}: duplicate sequence identity")
    if not np.isfinite(vectors).all():
        raise ValueError(f"{path}: non-finite vector")
    return dict(zip(ids, vectors))


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def hamming(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("RNA01-02 ladder sequences must have equal length")
    return sum(a != b for a, b in zip(left, right))


def validate_refset(
    object_rows: list[dict], edit_rows: list[dict], population: dict
) -> dict:
    objects: dict[str, dict] = {}
    for row in object_rows:
        object_id = row.get("object_id")
        sequence = row.get("sequence")
        sequence_identity = row.get("sequence_sha256")
        if not isinstance(object_id, str) or object_id in objects:
            raise ValueError("duplicate or invalid object_id in OBJECT_LEDGER.json")
        if not isinstance(sequence, str) or not sequence or not set(sequence).issubset(CANONICAL_RNA):
            raise ValueError(f"{object_id}: invalid canonical RNA sequence")
        if row.get("sequence_length_nt") != len(sequence):
            raise ValueError(f"{object_id}: sequence length metadata mismatch")
        if sequence_sha256(sequence) != sequence_identity:
            raise ValueError(f"{object_id}: sequence SHA256 mismatch")
        objects[object_id] = row

    leaves: dict[str, dict] = {}
    referenced_objects: set[str] = set()
    for row in edit_rows:
        leaf_id = row.get("leaf_id")
        if not isinstance(leaf_id, str) or leaf_id in leaves:
            raise ValueError("duplicate or invalid leaf_id in EDIT_LADDER_LEDGER.json")
        resolved: dict[str, dict] = {}
        for state in ("x0", "x1", "x2"):
            object_id = row.get(f"{state}_object_id")
            if object_id not in objects:
                raise ValueError(f"{leaf_id}: absent {state} object")
            obj = objects[object_id]
            if obj["sequence_sha256"] != row.get(f"{state}_sequence_sha256"):
                raise ValueError(f"{leaf_id}: {state} object/ladder sequence identity mismatch")
            resolved[state] = obj
            referenced_objects.add(object_id)
        x0, x1, x2 = (resolved[state]["sequence"] for state in ("x0", "x1", "x2"))
        distances = (hamming(x0, x1), hamming(x1, x2), hamming(x0, x2))
        if distances != (1, 1, 2):
            raise ValueError(f"{leaf_id}: expected exact Hamming 1/1/2 ladder")
        if tuple(row.get(key) for key in ("hamming_01", "hamming_12", "hamming_02")) != distances:
            raise ValueError(f"{leaf_id}: declared Hamming metadata mismatch")
        a = row.get("edit_a_position_0based")
        b = row.get("edit_b_position_0based")
        if not isinstance(a, int) or not isinstance(b, int) or a == b or not (0 <= a < len(x0)) or not (0 <= b < len(x0)):
            raise ValueError(f"{leaf_id}: invalid edit coordinates")
        if (x0[a], x1[a]) != (row.get("edit_a_from"), row.get("edit_a_to")):
            raise ValueError(f"{leaf_id}: first-edit identity mismatch")
        if (x1[b], x2[b]) != (row.get("edit_b_from"), row.get("edit_b_to")):
            raise ValueError(f"{leaf_id}: second-edit identity mismatch")
        leaves[leaf_id] = row

    if set(objects) != referenced_objects:
        raise ValueError("OBJECT_LEDGER.json contains missing or extra formal objects")
    members = population.get("members")
    if not isinstance(members, list) or len(members) != len(edit_rows):
        raise ValueError("population member count mismatch")
    projected = [{key: leaves[row["leaf_id"]].get(key) for key in POPULATION_MEMBER_KEYS} for row in members]
    if any(set(row) != set(POPULATION_MEMBER_KEYS) for row in members):
        raise ValueError("population member schema mismatch")
    if projected != members:
        raise ValueError("population members do not match the edit ladder ledger")
    if len({row["leaf_id"] for row in members}) != len(members):
        raise ValueError("duplicate leaf_id in population manifest")
    actual_clusters = len({row["outer_cluster_id"] for row in members})
    if population.get("outer_cluster_count") != actual_clusters:
        raise ValueError("population outer-cluster count mismatch")
    canonical_members = sorted(members, key=lambda row: row["leaf_id"])
    identity = hashlib.sha256(
        json.dumps(canonical_members, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()
    declared_identity = population.get("identity_sha256")
    if declared_identity != identity:
        raise ValueError("population identity declaration does not match actual content")
    return {
        "identity_sha256": identity,
        "registered_population_identity_match": identity == REGISTERED_POPULATION_IDENTITY_SHA256,
        "outer_cluster_count": actual_clusters,
    }


def safe_distance(left: np.ndarray, right: np.ndarray) -> float:
    left_norm = float(np.linalg.norm(left))
    right_norm = float(np.linalg.norm(right))
    if left_norm == 0.0 or right_norm == 0.0:
        return 1.0
    return 1.0 - float(np.dot(left, right) / (left_norm * right_norm))


def contrasts(row: dict, vectors: dict[str, np.ndarray]) -> tuple[float, float, float, float, float]:
    x0 = vectors[row["x0_sequence_sha256"]]
    x1 = vectors[row["x1_sequence_sha256"]]
    x2 = vectors[row["x2_sequence_sha256"]]
    d01 = safe_distance(x0, x1)
    d12 = safe_distance(x1, x2)
    d02 = safe_distance(x0, x2)
    return d01, d12, d02, d02 - d01, d02 - d12


def required_ids(rows: list[dict]) -> set[str]:
    return {row[key] for row in rows for key in ("x0_sequence_sha256", "x1_sequence_sha256", "x2_sequence_sha256")}


def require_exact_coverage(label: str, vectors: dict[str, np.ndarray], required: set[str]) -> None:
    actual = set(vectors)
    if actual != required:
        raise ValueError(f"{label}: ID coverage mismatch; missing={len(required - actual)}, extra={len(actual - required)}")


def derive_tau(rows: list[dict], first: dict[str, np.ndarray], second: dict[str, np.ndarray]) -> float:
    deviations: list[float] = []
    for row in rows:
        *_, left_one, right_one = contrasts(row, first)
        *_, left_two, right_two = contrasts(row, second)
        deviations.extend((abs(left_one - left_two), abs(right_one - right_two)))
    return float(max(deviations))


def credit(delta: float, tau: float) -> float:
    if delta > tau:
        return 1.0
    if delta < -tau:
        return -1.0
    return 0.0


def bootstrap(values: np.ndarray, seed: int) -> tuple[list[float], np.ndarray]:
    rng = np.random.default_rng(seed)
    count = len(values)
    means = np.empty(BOOTSTRAP_DRAWS, dtype=np.float64)
    for start in range(0, BOOTSTRAP_DRAWS, 2_000):
        size = min(2_000, BOOTSTRAP_DRAWS - start)
        means[start : start + size] = values[rng.integers(0, count, size=(size, count))].mean(axis=1)
    return [float(value) for value in np.quantile(means, [0.025, 0.975])], means


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="RNA01-02 L16K public scorer")
    parser.add_argument("--refset-dir", type=Path, default=Path(__file__).resolve().parents[1] / "refset")
    parser.add_argument("--vectors", type=Path, required=True)
    parser.add_argument("--calibration-repeat-1", type=Path, required=True)
    parser.add_argument("--calibration-repeat-2", type=Path, required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    edit_root = args.refset_dir / "RNA0102_L16K_EDIT_REFSET"
    object_path = edit_root / "OBJECT_LEDGER.json"
    edit_path = edit_root / "EDIT_LADDER_LEDGER.json"
    population_path = edit_root / "PLANNED_POPULATION_MANIFEST.json"
    object_rows = load_json(object_path)
    edit_rows = load_json(edit_path)
    calibration_rows = load_json(args.refset_dir / "RNA0102_calibration_ladders.json")
    population = load_json(population_path)
    refset_identity = validate_refset(object_rows, edit_rows, population)
    formal_by_leaf = {row["leaf_id"]: row for row in edit_rows}
    if len(calibration_rows) != 128 or len({row.get("leaf_id") for row in calibration_rows}) != 128:
        raise ValueError("calibration panel must contain 128 unique formal ladders")
    if any(formal_by_leaf.get(row.get("leaf_id")) != row for row in calibration_rows):
        raise ValueError("calibration panel does not match the validated formal ladder ledger")
    vectors = load_vectors(args.vectors)
    repeat_one = load_vectors(args.calibration_repeat_1)
    repeat_two = load_vectors(args.calibration_repeat_2)
    require_exact_coverage("formal vectors", vectors, required_ids(edit_rows))
    calibration_ids = required_ids(calibration_rows)
    require_exact_coverage("calibration repeat 1", repeat_one, calibration_ids)
    require_exact_coverage("calibration repeat 2", repeat_two, calibration_ids)
    formal_dimension = next(iter(vectors.values())).shape
    if next(iter(repeat_one.values())).shape != formal_dimension or next(iter(repeat_two.values())).shape != formal_dimension:
        raise ValueError("formal and calibration representation dimensions differ")
    tau = derive_tau(calibration_rows, repeat_one, repeat_two)

    individual_rows: list[dict] = []
    by_cluster: dict[str, list[float]] = {}
    for row in edit_rows:
        d01, d12, d02, left_delta, right_delta = contrasts(row, vectors)
        left_credit = credit(left_delta, tau)
        right_credit = credit(right_delta, tau)
        score = (left_credit + right_credit) / 2.0
        cluster_id = row["outer_cluster_id"]
        by_cluster.setdefault(cluster_id, []).append(score)
        individual_rows.append({
            "leaf_id": row["leaf_id"], "outer_cluster_id": cluster_id,
            "x0_object_id": row["x0_object_id"], "x1_object_id": row["x1_object_id"], "x2_object_id": row["x2_object_id"],
            "d01": format(d01, ".17g"), "d12": format(d12, ".17g"), "d02": format(d02, ".17g"),
            "left_delta": format(left_delta, ".17g"), "right_delta": format(right_delta, ".17g"),
            "left_F_credit": format(left_credit, ".1f"), "right_F_credit": format(right_credit, ".1f"),
            "individual_F": format(score, ".2f"),
        })

    cluster_rows = [{
        "outer_cluster_id": cluster_id,
        "leaf_count": len(by_cluster[cluster_id]),
        "cluster_F": format(float(np.mean(by_cluster[cluster_id], dtype=np.float64)), ".17g"),
    } for cluster_id in sorted(by_cluster)]
    cluster_values = np.asarray([float(row["cluster_F"]) for row in cluster_rows], dtype=np.float64)
    if len(cluster_rows) != refset_identity["outer_cluster_count"]:
        raise ValueError("outer-cluster count mismatch")
    seed = int(hashlib.sha256(f"RNA0102|{args.model_id}|L16K|outer-cluster-bootstrap-v1".encode()).hexdigest()[:16], 16)
    bounds, bootstrap_means = bootstrap(cluster_values, seed)
    formal_F = float(cluster_values.mean())

    args.output.mkdir(parents=True, exist_ok=True)
    individual_path = args.output / "INDIVIDUAL_SCORES.tsv"
    cluster_path = args.output / "CLUSTER_SCORE_LEDGER.tsv"
    calibration_path = args.output / "CALIBRATION.json"
    result_path = args.output / "FORMAL_RESULT.json"
    write_tsv(individual_path, list(individual_rows[0]), individual_rows)
    write_tsv(cluster_path, list(cluster_rows[0]), cluster_rows)
    calibration_path.write_text(json.dumps({
        "status": "PASS", "model_id": args.model_id, "calibration_ladders": len(calibration_rows),
        "complete_shared_endpoint_contrast_deviations": len(calibration_rows) * 2,
        "tau_delta": tau, "policy": CALIBRATION_POLICY, "deterministic": tau == 0.0,
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    result_path.write_text(json.dumps({
        "status": "PASS", "formal_result": True, "assay_id": "C-RNA-01-02", "model_id": args.model_id,
        "population": "L16K" if refset_identity["registered_population_identity_match"] else population.get("population_id", "UNREGISTERED"),
        "reference_population_identity_sha256": refset_identity["identity_sha256"],
        "registered_reference_population": refset_identity["registered_population_identity_match"],
        "evaluation_leaves": len(individual_rows), "outer_clusters": len(cluster_rows), "tau_delta": tau,
        "formal_interface": "F_ONLY", "score_transform_version": "RNA01-02_v0.2.6_F_ONLY_DIRECT_CREDITS",
        "metric": "F_RNA01-02", "F": formal_F, "F_95_CI": bounds,
        "bootstrap_draws": BOOTSTRAP_DRAWS, "bootstrap_unit": "outer_cluster", "bootstrap_seed": seed,
        "coverage": {"required": len(required_ids(edit_rows)), "scored": len(vectors)},
        "input_identities": {
            "object_ledger": sha256(object_path), "edit_ladder_ledger": sha256(edit_path),
            "population_manifest": sha256(population_path), "vectors": sha256(args.vectors),
            "calibration_repeat_1": sha256(args.calibration_repeat_1),
            "calibration_repeat_2": sha256(args.calibration_repeat_2),
        },
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    np.save(args.output / "BOOTSTRAP_MEANS.npy", bootstrap_means)
    outputs = [individual_path, cluster_path, calibration_path, result_path, args.output / "BOOTSTRAP_MEANS.npy"]
    (args.output / "OUTPUT_MANIFEST.json").write_text(json.dumps({
        "files": [{"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs]
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
