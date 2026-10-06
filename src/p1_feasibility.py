"""P1 Phase-2 feasibility statistics + dataset-version export.

Run: python3 src/p1_feasibility.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import STANDARDIZED, SPLITS, RESULTS, PROJECTS

P1 = PROJECTS / "P1_rr_benchmark"


def load():
    return pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")


def main():
    df = load()
    usable = df[df.ln_r1.notna() & df.ln_r2.notna()].copy()
    prim = usable[usable.is_primary_measurement].copy()

    # --- dataset versions -------------------------------------------------
    version_u = prim.copy()
    version_hq = prim[prim.quality_grade.isin(["A", "B"])].copy()
    version_c = prim[(prim.temperature_K.notna())
                     & (prim.solvent_family != "unknown")].copy()
    version_legacy = usable.copy()

    for name, v in [("dataset_U_unique", version_u),
                    ("dataset_HQ_highquality", version_hq),
                    ("dataset_C_condition", version_c),
                    ("dataset_legacy", version_legacy)]:
        v.to_parquet(SPLITS / f"{name}.parquet", index=False)

    stats = {
        "N_raw_measurements": int(len(df)),
        "N_unique_measurements_dedup": int(len(version_u)),
        "N_unique_pairs": int(usable.canonical_pair_id.nunique()),
        "N_unique_monomers": int(pd.unique(
            pd.concat([usable.m1_inchikey, usable.m2_inchikey]).dropna()).size),
        "N_condition_complete": int(len(version_c)),
        "N_quality_A_B": int(len(version_hq)),
        "N_quality_A": int((prim.quality_grade == "A").sum()),
        "N_quality_B": int((prim.quality_grade == "B").sum()),
        "N_quality_C": int((prim.quality_grade == "C").sum()),
        "N_quality_D": int((prim.quality_grade == "D").sum()),
        "N_dup_groups_multi_source":
            int((df.groupby("duplicate_group_id").source_dataset
                 .nunique() > 1).sum()),
        "pct_condition_complete":
            round(100 * len(version_c) / max(len(version_u), 1), 1),
    }

    # --- P2 multi-temperature feasibility ---------------------------------
    sub = usable[usable.temperature_K.notna()].copy()
    sub["pair_solvent"] = (sub.canonical_pair_id + "|" + sub.solvent_family)
    grp = sub.groupby("pair_solvent")
    n_temp = grp.temperature_K.nunique()
    span = grp.temperature_K.apply(lambda s: s.max() - s.min())
    stats["P2_groups_ge3_temperatures"] = int((n_temp >= 3).sum())
    stats["P2_groups_ge4_temperatures"] = int((n_temp >= 4).sum())
    stats["P2_median_temperature_span_K"] = (
        round(float(span[n_temp >= 3].median()), 1)
        if (n_temp >= 3).any() else None)
    stats["P2_groups_ge3_with_span_ge30K"] = int(
        ((n_temp >= 3) & (span >= 30)).sum())

    print(json.dumps(stats, indent=2))
    (P1 / "feasibility_stats.json").write_text(json.dumps(stats, indent=2))

    # Go/no-go assessment
    go = (stats["N_unique_measurements_dedup"] >= 1500
          and stats["N_condition_complete"] >= 800)
    print("\nP1 Go/No-Go (>=1500 unique & >=800 condition-complete):",
          "GO" if go else "CONDITIONAL")
    p2_go = stats["P2_groups_ge3_temperatures"]
    print(f"P2 strict multi-T feasibility: {p2_go} groups (>=150 strong Go, "
          f"50-150 conditional, <50 no-go)")


if __name__ == "__main__":
    main()
