#!/usr/bin/env python3
"""Score the registered RNApaln measurements through the RNA09 core."""
import argparse, csv, hashlib, importlib.util, json
from pathlib import Path
import numpy as np


def load_public_core():
    path = Path(__file__).resolve().parents[2] / "scoring/score_official.py"
    spec = importlib.util.spec_from_file_location("rna09_public_scorer", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


parser = argparse.ArgumentParser(); parser.add_argument("--panel", type=Path, required=True); parser.add_argument("--score-dir", type=Path, required=True)
args = parser.parse_args(); args.score_dir.mkdir(parents=True, exist_ok=True); core = load_public_core()
scores = json.loads((args.score_dir / "RAW_DIRECTIONAL_SCORES.json").read_text())
with (args.panel / "scoring_input_v01/RNA09_ENDPOINT_OBJECT_MAP.tsv").open() as handle:
    object_map = {row["object_id"]: row["sequence_sha256"] for row in csv.DictReader(handle, delimiter="\t")}
rows = core.read_tsv(args.panel / "FORMAL_CONTRAST_LEDGER.tsv")


def similarity(a, b):
    raw = lambda left, right: scores[left + "\t" + right]
    denominator = raw(a, a) + raw(b, b)
    if denominator <= 0: raise ValueError("nonpositive self denominator")
    return (raw(a, b) + raw(b, a)) / denominator


for row in rows:
    anchor, positive, negative = [object_map[row[key]] for key in ("anchor_id", "positive_partner_id", "negative_partner_id")]
    row["positive_similarity"], row["negative_similarity"] = similarity(anchor, positive), similarity(anchor, negative)
core.assign_credits(rows); aggregation = core.aggregate_relation_rows(rows)
leaf_path = args.score_dir / "RNAPALN_PRIMARY_LEAF_LEDGER.tsv"; core.write_tsv(leaf_path, rows)
core.write_hierarchy_outputs(args.score_dir, aggregation, "RNAPALN_PRIMARY")
np.save(args.score_dir / "RNAPALN_PRIMARY_F_BOOTSTRAP_MEANS.npy", aggregation["F_bootstrap"])
result = {"status":"PASS", "contrasts":len(rows), "anchor_bundles":len(aggregation["anchor_scores"]),
          "source_parents":len(aggregation["parent_scores"]), "source_dependency_units":len(aggregation["unit_scores"]),
          "structure_relation_groups":1, "F":aggregation["F"], "F_95_CI":aggregation["F_95_CI"],
          "formal_interface":"F_ONLY", "score_transform_version":"RNA09_F_ONLY_v1", "bootstrap_seed":core.BOOTSTRAP_SEED, "bootstrap_draws":core.BOOTSTRAP_DRAWS,
          "leaf_ledger_sha256":hashlib.sha256(leaf_path.read_bytes()).hexdigest(),
          "raw_score_sha256":hashlib.sha256((args.score_dir / "RAW_DIRECTIONAL_SCORES.json").read_bytes()).hexdigest()}
(args.score_dir / "RNAPALN_PRIMARY_RESULT.json").write_text(json.dumps(result, indent=2)+"\n"); print(json.dumps(result))
