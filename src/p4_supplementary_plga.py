"""Composition-matched sequence/Tg check on published PLGA measurements.

Source: Patil & Won, Macromolecules 2024, DOI 10.1021/acs.macromol.4c00106.
The extracted 12-row input is stored beside the project outputs so this check
does not depend on the sibling SequenceTg repository.

This is an exploratory within-composition association analysis, not a causal
or cross-chemistry validation. It uses exact permutations of Tg labels within
each of the four composition groups.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projects" / "P4_sequence_tg"
INPUT = PROJECT / "plga_sequence_tg_input.csv"


def within_group_center(df: pd.DataFrame, column: str) -> np.ndarray:
    return (df[column] - df.groupby("wL")[column].transform("mean")).to_numpy(float)


def exact_group_permutation(df: pd.DataFrame, feature: str) -> dict:
    """Two-sided exact permutation test, permuting Tg only within composition."""
    x = within_group_center(df, feature)
    y = within_group_center(df, "Tg_C")
    observed = float(np.corrcoef(x, y)[0, 1])
    groups = [g.index.to_list() for _, g in df.groupby("wL", sort=True)]
    group_perms = [list(itertools.permutations(idx)) for idx in groups]
    y_values = df.Tg_C.to_numpy(float)
    null = []
    for choice in itertools.product(*group_perms):
        yp = np.empty(len(df), dtype=float)
        for indices, permuted_indices in zip(groups, choice):
            for target, source in zip(indices, permuted_indices):
                yp[target] = y_values[source]
        for indices in groups:
            yp[indices] -= yp[indices].mean()
        null.append(float(np.corrcoef(x, yp)[0, 1]))
    null = np.asarray(null)
    slope = float(np.dot(x, y) / np.dot(x, x))
    return {
        "feature": feature,
        "within_composition_pearson_r": observed,
        "within_composition_slope_C_per_feature_unit": slope,
        "exact_two_sided_p": float(np.mean(np.abs(null) >= abs(observed))),
        "permutations": int(len(null)),
    }


def main() -> None:
    df = pd.read_csv(INPUT).sort_values(["wL", "feed_rate"], ascending=[True, False])
    features = ["blockiness", "alternation", "LL", "LG", "feed_rate", "dispersity"]
    tests = [exact_group_permutation(df.reset_index(drop=True), c) for c in features]
    # The sequence, process, and dispersity associations are exploratory as a
    # family; report Holm-adjusted p-values to discourage selective emphasis.
    order = np.argsort([r["exact_two_sided_p"] for r in tests])
    adjusted = np.empty(len(tests), dtype=float)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (len(tests) - rank) * tests[idx]["exact_two_sided_p"])
        adjusted[idx] = min(1.0, running)
    for row, p_adj in zip(tests, adjusted):
        row["holm_adjusted_p"] = float(p_adj)

    contrasts = []
    for wL, group in df.groupby("wL", sort=True):
        high = group.loc[group.feed_rate.idxmax()]
        low = group.loc[group.feed_rate.idxmin()]
        contrasts.append({
            "wL": float(wL),
            "high_feed_rate": float(high.feed_rate),
            "low_feed_rate": float(low.feed_rate),
            "delta_Tg_low_minus_high_C": float(low.Tg_C - high.Tg_C),
            "delta_blockiness_low_minus_high": float(low.blockiness - high.blockiness),
            "delta_alternation_low_minus_high": float(low.alternation - high.alternation),
            "delta_dispersity_low_minus_high": float(low.dispersity - high.dispersity),
        })

    out = {
        "source": "Patil & Won, Macromolecules 2024, DOI 10.1021/acs.macromol.4c00106",
        "n_samples": int(len(df)),
        "n_composition_groups": int(df.wL.nunique()),
        "analysis": "exact within-composition permutation; exploratory, one chemistry pair",
        "tests": tests,
        "matched_extreme_feed_contrasts": contrasts,
        "limitations": [
            "one monomer pair and four composition groups only",
            "feed rate changes sequence and may change other polymerization outcomes",
            "molecular-weight dispersity is available, but absolute molecular weights are not in this extracted table",
            "sequence/Tg observations are from one published study, not a prospective validation",
        ],
    }
    (PROJECT / "p4_supplementary_plga.json").write_text(json.dumps(out, indent=2))
    pd.DataFrame(tests).to_csv(PROJECT / "p4_supplementary_plga_tests.csv", index=False)
    pd.DataFrame(contrasts).to_csv(PROJECT / "p4_supplementary_plga_contrasts.csv", index=False)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
