#!/usr/bin/env python3
"""Fixed-panel test of knowledge-specific unevenness.

The estimand is the mean formal F across the five prespecified models for each
Constraint. The null states that, within a biological domain, Constraint labels
do not identify a shared high- or low-accessibility profile across models.

Calibration independently permutes Constraint labels within each model and
domain. This preserves every model's observed F distribution and the domain-
specific model panels while destroying alignment on a common Constraint label.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "Supplementary_Table_S4_Complete_Model_Performance.tsv"
OUT = ROOT / "results"

DOMAIN_MAP = {
    "DNA01": "DNA", "DNA02-02": "DNA", "DNA03": "DNA", "DNA06": "DNA",
    "RNA01-02": "RNA", "RNA04": "RNA", "RNA05": "RNA", "RNA06": "RNA",
    "RNA09": "RNA", "RNA10": "RNA",
    "PS01": "Protein sequence", "PS02A": "Protein sequence", "PS02B": "Protein sequence",
    "PS04": "Protein sequence", "PS05": "Protein sequence", "PS08": "Protein sequence",
    "PS09": "Protein sequence", "PS10": "Protein sequence",
    "ST01": "Protein structure", "ST02": "Protein structure", "ST04": "Protein structure",
    "ST05": "Protein structure",
}
DOMAIN_ORDER = ["DNA", "RNA", "Protein sequence", "Protein structure"]
SEED = 20261025
REPLICATES = 1_000_000
BATCH_SIZE = 2_000


def load_panel() -> list[dict[str, str]]:
    with INPUT.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 110:
        raise ValueError(f"Expected 110 model-by-Constraint cells, found {len(rows)}")
    for row in rows:
        row["Domain"] = DOMAIN_MAP.get(row["Constraint"], "")
        if not row["Domain"]:
            raise ValueError(f"Unmapped Constraint: {row['Constraint']}")
        row["F_value"] = float(row["F"])
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["Constraint"]] = counts.get(row["Constraint"], 0) + 1
    if set(counts.values()) != {5}:
        raise ValueError(f"Each Constraint must have five models: {counts}")
    return rows


def domain_matrices(rows: list[dict[str, str]]) -> dict[str, tuple[list[str], list[str], np.ndarray]]:
    result = {}
    for domain in DOMAIN_ORDER:
        subset = [r for r in rows if r["Domain"] == domain]
        constraints = list(dict.fromkeys(r["Constraint"] for r in subset))
        models = list(dict.fromkeys(r["Model"] for r in subset))
        matrix = np.empty((len(constraints), len(models)), dtype=np.float64)
        lookup = {(r["Constraint"], r["Model"]): r["F_value"] for r in subset}
        for i, constraint in enumerate(constraints):
            for j, model in enumerate(models):
                matrix[i, j] = lookup[(constraint, model)]
        if matrix.shape[1] != 5:
            raise ValueError(f"{domain} does not contain five models: {matrix.shape}")
        result[domain] = (constraints, models, matrix)
    return result


def knowledge_statistic(matrix: np.ndarray) -> float:
    constraint_means = matrix.mean(axis=1)
    return float(np.sum((constraint_means - constraint_means.mean()) ** 2))


def permuted_statistics(matrix: np.ndarray, draws: int, rng: np.random.Generator) -> np.ndarray:
    constraints, models = matrix.shape
    random_keys = rng.random((draws, models, constraints))
    permutations = np.argsort(random_keys, axis=2)
    permuted = np.empty((draws, constraints, models), dtype=np.float64)
    for model_index in range(models):
        permuted[:, :, model_index] = matrix[permutations[:, model_index, :], model_index]
    means = permuted.mean(axis=2)
    centered = means - means.mean(axis=1, keepdims=True)
    return np.sum(centered * centered, axis=1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = load_panel()
    matrices = domain_matrices(rows)

    observed_by_domain = []
    constraint_summary = []
    observed_total = 0.0
    for domain in DOMAIN_ORDER:
        constraints, _models, matrix = matrices[domain]
        statistic = knowledge_statistic(matrix)
        observed_total += statistic
        observed_by_domain.append((domain, len(constraints), statistic))
        for constraint, values in zip(constraints, matrix, strict=True):
            matching = next(r for r in rows if r["Constraint"] == constraint)
            constraint_summary.append({
                "Domain": domain,
                "Constraint": constraint,
                "Short_name": matching["Short_name"],
                "mean_F_across_models": float(values.mean()),
                "sd_F_across_models": float(values.std(ddof=1)),
                "min_F": float(values.min()),
                "max_F": float(values.max()),
            })

    rng = np.random.default_rng(SEED)
    null = np.empty(REPLICATES, dtype=np.float64)
    for start in range(0, REPLICATES, BATCH_SIZE):
        stop = min(start + BATCH_SIZE, REPLICATES)
        batch = np.zeros(stop - start, dtype=np.float64)
        for domain in DOMAIN_ORDER:
            batch += permuted_statistics(matrices[domain][2], stop - start, rng)
        null[start:stop] = batch

    exceedances = int(np.count_nonzero(null >= observed_total))
    p_value = (exceedances + 1) / (REPLICATES + 1)
    q025, median, q975 = np.quantile(null, [0.025, 0.5, 0.975])

    with (OUT / "constraint_mean_F_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(constraint_summary[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(constraint_summary)

    with (OUT / "knowledge_unevenness_by_domain.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["domain", "constraints", "observed_constraint_mean_SS"])
        writer.writerows(observed_by_domain)

    with (OUT / "global_knowledge_unevenness_permutation.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "constraints", "formal_F_cells", "domains", "permutation_replicates", "seed",
            "observed_constraint_mean_SS", "null_median", "null_q025", "null_q975",
            "exceedances", "p_value",
        ])
        writer.writerow([
            22, len(rows), 4, REPLICATES, SEED, observed_total, median, q025, q975,
            exceedances, p_value,
        ])

    np.savetxt(
        OUT / "global_knowledge_unevenness_null_draws.tsv",
        null,
        delimiter="\t",
        header="null_constraint_mean_SS",
        comments="",
        fmt="%.12g",
    )

    p_text = f"P = {p_value:.6g}" if exceedances else f"P < {1 / (REPLICATES + 1):.2g}"
    summary = f"""# Global knowledge-specific unevenness test

## Estimand

For each of the 22 prespecified Constraints, accessibility is summarized as the mean formal F across that domain's fixed panel of five models. The test asks whether these Constraint means are more uneven than expected if Constraint identity carried no shared high- or low-accessibility signal across models.

## Null hypothesis

Within each biological domain, Constraint labels are exchangeable independently for each model. Each model's observed F distribution and each domain's five-model panel are fixed, but alignment of scores on a common biological Constraint is removed.

## Statistic and calibration

The statistic is the sum, across domains and Constraints, of squared deviations of each five-model mean F from its domain mean. Each of the 22 Constraints therefore contributes once. Calibration uses {REPLICATES:,} independent within-model, within-domain label permutations (seed {SEED}).

Observed knowledge unevenness SS = {observed_total:.8f}; null median = {median:.8f}; 95% null interval = [{q025:.8f}, {q975:.8f}]; exceedances = {exceedances}; **{p_text}**.

## Interpretation boundary

This is a knowledge/Constraint main-effect test for the fixed model panels: it tests whether some measured biological relations are consistently more accessible than others when averaged across the five models in their domain. It is distinct from the model-by-Constraint interaction test, which asks whether relative model advantage changes across Constraints. It does not infer an unobserved superpopulation of all possible biological knowledge and does not identify intrinsic biological difficulty independently of the chosen models, RefSets and measurement protocols.

No Forward pass, primary Scorer run or formal F value was changed.
"""
    (OUT / "GLOBAL_KNOWLEDGE_UNEVENNESS_SUMMARY.md").write_text(summary, encoding="utf-8")
    print(
        f"observed={observed_total:.10f} null_median={median:.10f} "
        f"exceedances={exceedances} p={p_value:.10g}"
    )


if __name__ == "__main__":
    main()
