from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


LAMBDA_GRID = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0)
BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 20261003


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("endpoint ledger is empty")
    return rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def training_weights(rows: list[dict[str, str]]) -> np.ndarray:
    hierarchy = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    for index, row in enumerate(rows):
        hierarchy[row["outer_cluster_id"]][row["parent_id"]][row["block_id"]][row["paired_state"]].append(index)
    result = np.zeros(len(rows))
    for parents in hierarchy.values():
        for blocks in parents.values():
            for states in blocks.values():
                for indices in states.values():
                    result[indices] = 1 / len(hierarchy) / len(parents) / len(blocks) / len(states) / len(indices)
    return result / result.mean()


def block_scores(rows: list[dict[str, str]], predictions: np.ndarray) -> tuple[dict[str, float], list[dict]]:
    grouped = defaultdict(lambda: defaultdict(list))
    identity: dict[str, tuple[str, str]] = {}
    for row, prediction in zip(rows, predictions):
        grouped[row["block_id"]][row["paired_state"]].append(float(prediction))
        current = (row["parent_id"], row["outer_cluster_id"])
        if row["block_id"] in identity and identity[row["block_id"]] != current:
            raise ValueError("block maps to multiple parent or cluster identities")
        identity[row["block_id"]] = current
    values: dict[str, float] = {}
    records: list[dict] = []
    for block_id in sorted(grouped):
        states = grouped[block_id]
        if not states["paired"] or not states["unpaired"]:
            raise ValueError(f"incomplete block {block_id}")
        paired_mean = float(np.mean(states["paired"]))
        unpaired_mean = float(np.mean(states["unpaired"]))
        score = 1.0 if paired_mean > unpaired_mean else -1.0 if paired_mean < unpaired_mean else 0.0
        values[block_id] = score
        parent, cluster = identity[block_id]
        records.append({
            "block_id": block_id, "parent_id": parent, "outer_cluster_id": cluster,
            "paired_count": len(states["paired"]), "unpaired_count": len(states["unpaired"]),
            "paired_mean_prediction": format(paired_mean, ".17g"),
            "unpaired_mean_prediction": format(unpaired_mean, ".17g"),
            "raw_contrast": format(paired_mean - unpaired_mean, ".17g"),
            "F_credit": format(score, ".17g"),
        })
    return values, records


def aggregate(rows: list[dict[str, str]], values: dict[str, float]) -> tuple[float, dict[str, float], dict[str, float]]:
    parent_blocks: dict[str, list[float]] = defaultdict(list)
    parent_cluster: dict[str, str] = {}
    seen: set[str] = set()
    for row in rows:
        block = row["block_id"]
        if block in seen:
            continue
        seen.add(block)
        parent = row["parent_id"]
        parent_blocks[parent].append(values[block])
        parent_cluster[parent] = row["outer_cluster_id"]
    parents = {key: float(np.mean(value)) for key, value in parent_blocks.items()}
    cluster_parents: dict[str, list[float]] = defaultdict(list)
    for parent, value in parents.items():
        cluster_parents[parent_cluster[parent]].append(value)
    clusters = {key: float(np.mean(value)) for key, value in cluster_parents.items()}
    return float(np.mean(list(clusters.values()))), parents, clusters


def deterministic_folds(groups: np.ndarray, count: int = 4) -> list[set[str]]:
    ordered = sorted(set(groups), key=lambda value: hashlib.sha256(value.encode()).hexdigest())
    return [{value for index, value in enumerate(ordered) if index % count == fold} for fold in range(count)]


def fit(train: np.ndarray, test: np.ndarray, features: np.ndarray, rows: list[dict[str, str]], penalty: float) -> np.ndarray:
    weights = training_weights([rows[index] for index in train])
    weights = weights / weights.sum()
    mean = (features[train] * weights[:, None]).sum(axis=0)
    scale = np.sqrt(((features[train] - mean) ** 2 * weights[:, None]).sum(axis=0))
    scale[scale == 0] = 1.0
    train_x = (features[train] - mean) / scale
    test_x = (features[test] - mean) / scale
    labels = np.asarray([rows[index]["paired_state"] == "paired" for index in train], dtype=float)

    def objective(theta: np.ndarray):
        linear = train_x @ theta[:-1] + theta[-1]
        loss = np.sum(weights * (np.logaddexp(0, linear) - labels * linear)) + 0.5 * penalty * np.dot(theta[:-1], theta[:-1])
        probability = 1 / (1 + np.exp(-np.clip(linear, -40, 40)))
        derivative = (probability - labels) * weights
        gradient = np.r_[train_x.T @ derivative + penalty * theta[:-1], derivative.sum()]
        return loss, gradient

    fitted = minimize(lambda theta: objective(theta), np.zeros(train_x.shape[1] + 1), jac=True, method="L-BFGS-B", options={"maxiter": 2000})
    if not fitted.success or not np.isfinite(fitted.x).all():
        raise RuntimeError(f"ridge optimizer failed: {fitted.message}")
    linear = test_x @ fitted.x[:-1] + fitted.x[-1]
    return 1 / (1 + np.exp(-np.clip(linear, -40, 40)))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records[current] = "".join(chunks).upper()
            current = line[1:].split()[0]
            if not current or current in records:
                raise ValueError("FASTA IDs must be nonempty and unique")
            chunks = []
        elif current is None:
            raise ValueError("sequence before first FASTA header")
        else:
            chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks).upper()
    if not records:
        raise ValueError("FASTA is empty")
    return records


def parse_position(value: str) -> int:
    if not isinstance(value, str) or re.fullmatch(r"0|[1-9][0-9]*", value) is None:
        raise ValueError(f"position0 must be a nonnegative integer literal, got {value!r}")
    return int(value)


def validate_refset(rows: list[dict[str, str]], blocks: list[dict[str, str]], fasta_path: Path) -> tuple[dict[str, str], str]:
    sequences = read_fasta(fasta_path)
    for sequence_id, sequence in sequences.items():
        if not sequence or set(sequence) - set("ACGU"):
            raise ValueError(f"invalid RNA sequence: {sequence_id}")
        if hashlib.sha256(sequence.encode()).hexdigest() != sequence_id:
            raise ValueError(f"FASTA sequence identity mismatch: {sequence_id}")
    block_table = {row["block_id"]: row for row in blocks}
    if len(block_table) != len(blocks):
        raise ValueError("duplicate block identity")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen_positions: set[tuple[str, int]] = set()
    parent_clusters: dict[str, str] = {}
    semantic_rows: list[dict[str, str]] = []
    for index, row in enumerate(rows):
        sequence_id = row["sequence_sha256"]
        if sequence_id not in sequences:
            raise ValueError(f"endpoint row {index} references a sequence outside FASTA")
        position = parse_position(row["position0"])
        sequence = sequences[sequence_id]
        if position >= len(sequence) or row["nucleotide_identity"] != sequence[position]:
            raise ValueError(f"endpoint biological coordinate mismatch at row {index}")
        key = (sequence_id, position)
        if key in seen_positions:
            raise ValueError("duplicate sequence-position endpoint")
        seen_positions.add(key)
        if row["paired_state"] not in {"paired", "unpaired"}:
            raise ValueError(f"invalid paired_state at row {index}")
        if row["relative_position_bin_count"] != "8" or not 0 <= parse_position(row["relative_position_bin"]) < 8:
            raise ValueError(f"invalid frozen 8-bin identity at row {index}")
        if row["complete_execution_bundle"].lower() != "true":
            raise ValueError(f"incomplete formal block at row {index}")
        parent, cluster = row["parent_id"], row["outer_cluster_id"]
        if parent in parent_clusters and parent_clusters[parent] != cluster:
            raise ValueError("parent maps to multiple dependency clusters")
        parent_clusters[parent] = cluster
        grouped[row["block_id"]].append(row)
        semantic_rows.append({key: row[key] for key in (
            "block_id", "sequence_sha256", "position0", "paired_state", "nucleotide_identity",
            "relative_position_bin", "parent_id", "outer_cluster_id",
        )})
    if set(grouped) != set(block_table):
        raise ValueError("endpoint and block ledgers have different block populations")
    identity_fields = (
        "record_id", "parent_id", "family_id", "outer_cluster_id", "encoded_input_hash",
        "primary_structure_map_id", "primary_condition_id", "nucleotide_identity",
        "relative_position_bin", "relative_position_bin_count",
    )
    for block_id, members in grouped.items():
        definition = block_table[block_id]
        for field in identity_fields:
            if len({row[field] for row in members}) != 1 or members[0][field] != definition[field]:
                raise ValueError(f"block membership mismatch for {block_id}: {field}")
        counts = {state: sum(row["paired_state"] == state for row in members) for state in ("paired", "unpaired")}
        if not all(counts.values()) or counts["paired"] != int(definition["paired_count"]) or counts["unpaired"] != int(definition["unpaired_count"]):
            raise ValueError(f"block state-count mismatch for {block_id}")
    payload = {"sequences": sorted(sequences), "endpoints": sorted(semantic_rows, key=lambda row: tuple(row.values()))}
    identity = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return sequences, identity


def load_features(path: Path, rows: list[dict[str, str]], released_sequences: dict[str, str]) -> tuple[np.ndarray, str]:
    archive = np.load(path, allow_pickle=False)
    required_positions = [(row["sequence_sha256"], parse_position(row["position0"])) for row in rows]
    if set(archive.files) == {"sequence_sha256", "position0", "vectors"}:
        raw_ids, raw_positions, raw_vectors = archive["sequence_sha256"], archive["position0"], archive["vectors"]
        if raw_ids.ndim != 1 or raw_ids.dtype.kind not in "US":
            raise ValueError("sequence_sha256 must be a one-dimensional string array")
        if raw_positions.ndim != 1 or raw_positions.dtype.kind not in "iu" or np.issubdtype(raw_positions.dtype, np.bool_):
            raise ValueError("position0 must be a one-dimensional integer array")
        if raw_vectors.dtype.kind not in "iuf" or np.issubdtype(raw_vectors.dtype, np.bool_):
            raise ValueError("vectors must contain real numeric values")
        if len(raw_ids) != len(raw_positions) or len(raw_ids) != len(raw_vectors):
            raise ValueError("sparse position arrays have different row counts")
        ids = [(sequence.decode() if isinstance(sequence, bytes) else str(sequence), int(position)) for sequence, position in zip(raw_ids, raw_positions)]
        vectors = np.asarray(raw_vectors, dtype=np.float64)
        if vectors.ndim != 2 or vectors.shape[1] <= 0 or len(ids) != len(vectors) or len(set(ids)) != len(ids) or not np.isfinite(vectors).all():
            raise ValueError("invalid sparse position-vector NPZ")
        if set(ids) != set(required_positions):
            raise ValueError(f"position coverage mismatch; missing={len(set(required_positions)-set(ids))}, extra={len(set(ids)-set(required_positions))}")
        lookup = {key: value for key, value in zip(ids, vectors)}
        return np.stack([lookup[key] for key in required_positions]), "sparse_position_npz"

    released_ids = set(released_sequences)
    required_sequences = {sequence for sequence, _ in required_positions}
    if set(archive.files) != released_ids:
        raise ValueError(f"sequence coverage mismatch; missing={len(released_ids-set(archive.files))}, extra={len(set(archive.files)-released_ids)}")
    if not required_sequences.issubset(released_ids):
        raise ValueError("endpoint ledger references sequences outside the released FASTA")
    picked = []
    dimension = None
    for sequence, position in required_positions:
        raw_matrix = archive[sequence]
        if raw_matrix.dtype.kind not in "iuf" or np.issubdtype(raw_matrix.dtype, np.bool_):
            raise ValueError(f"sequence matrix must contain real numeric values for {sequence}")
        matrix = np.asarray(raw_matrix, dtype=np.float64)
        if matrix.ndim != 2 or matrix.shape[0] != len(released_sequences[sequence]) or matrix.shape[1] <= 0 or not np.isfinite(matrix).all():
            raise ValueError(f"invalid sequence-position matrix for {sequence}")
        dimension = matrix.shape[1] if dimension is None else dimension
        if matrix.shape[1] != dimension:
            raise ValueError("representation dimension differs across sequences")
        picked.append(matrix[position])
    return np.stack(picked), "sequence_keyed_npz"


def main() -> None:
    parser = argparse.ArgumentParser(description="RNA05 public frozen nested-cross-fit scorer")
    default_ledger = Path(__file__).resolve().parents[1] / "refset" / "FORMAL_8BIN_ENDPOINT_LEDGER.tsv"
    parser.add_argument("--ledger", type=Path, default=default_ledger)
    refset = Path(__file__).resolve().parents[1] / "refset"
    default_fasta = refset / "RNA05_06_unique_inputs.fasta"
    default_blocks = refset / "RNA05_PRIMARY_MATCHED_BLOCKS_8BINS.tsv"
    parser.add_argument("--vectors", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, default=default_fasta)
    parser.add_argument("--blocks", type=Path, default=default_blocks)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = read_tsv(args.ledger)
    sequences, population_identity = validate_refset(rows, read_tsv(args.blocks), args.fasta)
    registered_identity = validate_refset(read_tsv(default_ledger), read_tsv(default_blocks), default_fasta)[1]
    features, representation_schema = load_features(args.vectors, rows, sequences)
    groups = np.asarray([row["outer_cluster_id"] for row in rows])
    predictions = np.full(len(rows), np.nan)
    selected_parameters = []
    for held_out in sorted(set(groups)):
        train = np.where(groups != held_out)[0]
        test = np.where(groups == held_out)[0]
        best = None
        for penalty in LAMBDA_GRID:
            validation_scores = []
            for inner_holdout in deterministic_folds(groups[train]):
                inner_train = train[np.asarray([group not in inner_holdout for group in groups[train]])]
                inner_test = train[np.asarray([group in inner_holdout for group in groups[train]])]
                inner_prediction = fit(inner_train, inner_test, features, rows, penalty)
                block_values, _ = block_scores([rows[index] for index in inner_test], inner_prediction)
                validation_scores.append(aggregate([rows[index] for index in inner_test], block_values)[0])
            key = (float(np.mean(validation_scores)), -penalty)
            if best is None or key > best[0]:
                best = (key, penalty)
        predictions[test] = fit(train, test, features, rows, best[1])
        selected_parameters.append({"held_out_cluster": held_out, "lambda": best[1], "inner_F": best[0][0]})
    if not np.isfinite(predictions).all():
        raise ValueError("missing held-out predictions")

    block_values, individual_rows = block_scores(rows, predictions)
    score, parent_values, cluster_values = aggregate(rows, block_values)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    # Preserve the frozen executed scorer's first-occurrence order.  Sorting
    # cluster IDs changes which score is associated with each seeded bootstrap
    # draw and therefore changes the otherwise identical percentile CI.
    cluster_ids = np.asarray(list(cluster_values))
    cluster_array = np.asarray([cluster_values[key] for key in cluster_ids])
    bootstrap = rng.choice(cluster_array, (BOOTSTRAP_DRAWS, len(cluster_array)), replace=True).mean(axis=1)
    low, high = np.quantile(bootstrap, [0.025, 0.975])

    args.output.mkdir(parents=True, exist_ok=True)
    endpoint_rows = [{
        "endpoint_index": index, "block_id": row["block_id"], "parent_id": row["parent_id"],
        "outer_cluster_id": row["outer_cluster_id"], "paired_state": row["paired_state"],
        "sequence_sha256": row["sequence_sha256"], "position0": row["position0"],
        "heldout_probability_paired": format(float(predictions[index]), ".17g"),
    } for index, row in enumerate(rows)]
    parent_rows = [{"parent_id": key, "parent_F": format(parent_values[key], ".17g")} for key in sorted(parent_values)]
    cluster_rows = [{"outer_cluster_id": key, "cluster_F": format(cluster_values[key], ".17g")} for key in sorted(cluster_values)]
    tables = {"HELDOUT_PREDICTIONS.tsv": endpoint_rows, "INDIVIDUAL_SCORES.tsv": individual_rows, "PARENT_SCORE_LEDGER.tsv": parent_rows, "CLUSTER_SCORE_LEDGER.tsv": cluster_rows}
    for name, records in tables.items():
        write_tsv(args.output / name, records)
    result_path = args.output / "FORMAL_RESULT.json"
    result_path.write_text(json.dumps({
        "status": "PASS", "formal_result": True, "assay_id": "C-RNA-05", "model_id": args.model_id,
        "representation_schema": representation_schema,
        "readout": "outer-cluster leave-one-out L2 ridge logistic regression with frozen inner selection",
        "endpoints": len(rows), "blocks": len(block_values), "parents": len(parent_values), "clusters": len(cluster_values),
        "F": score, "F_95_CI": [float(low), float(high)],
        "formal_interface": "RNA05_F_ONLY_v0.1.3", "score_transform_version": "direct_F_credit_v1",
        "outer_LOO": selected_parameters,
        "bootstrap_replicates": BOOTSTRAP_DRAWS, "bootstrap_seed": BOOTSTRAP_SEED, "bootstrap_unit": "outer_cluster_F",
        "ledger_sha256": sha256(args.ledger), "blocks_sha256": sha256(args.blocks), "fasta_sha256": sha256(args.fasta),
        "population_identity_sha256": population_identity, "registered_reference_population": population_identity == registered_identity,
        "vectors_sha256": sha256(args.vectors),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    np.save(args.output / "BOOTSTRAP_MEANS.npy", bootstrap)
    outputs = [args.output / name for name in tables] + [result_path, args.output / "BOOTSTRAP_MEANS.npy"]
    (args.output / "OUTPUT_MANIFEST.json").write_text(json.dumps({"files": [
        {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs
    ]}, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
