"""P1: gold-subset experiment — raise R² by aggregating repeated measurements.

Key idea: the R² ceiling is set by within-pair measurement noise (~38% of
variance, ceiling ~0.62). For pairs with >=3 independent measurements, the
pair-MEDIAN is a much lower-noise estimate of the true ln r, so the noise
floor rises and a model should reach higher R² at lower MAE.

Run: python3 src/p1_gold_subset.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import STANDARDIZED, RESULTS, PROJECTS
import p1_features as F

P1 = PROJECTS / "P1_rr_benchmark"


def load_primary():
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    df = df[df.is_primary_measurement &
            df.ln_r1.notna() & df.ln_r2.notna() &
            (df.r1 > 0) & (df.r2 > 0)].reset_index(drop=True)
    return df


def pair_medians(df):
    g = df.groupby("canonical_pair_id")
    return pd.DataFrame({
        "ln_r1": g.ln_r1.median(),
        "ln_r2": g.ln_r2.median(),
        "ln_r1_std": g.ln_r1.std(),
        "ln_r2_std": g.ln_r2.std(),
        "n_meas": g.size(),
        "m1_inchikey": g.m1_inchikey.first(),
        "m2_inchikey": g.m2_inchikey.first(),
    }).reset_index()


def struct_feats(pair_df, ee, dd):
    e1 = ee.reindex(pair_df.m1_inchikey).to_numpy(np.float32)
    e2 = ee.reindex(pair_df.m2_inchikey).to_numpy(np.float32)
    d1 = dd.reindex(pair_df.m1_inchikey).to_numpy(np.float32)
    d2 = dd.reindex(pair_df.m2_inchikey).to_numpy(np.float32)
    return np.hstack([e1, e2, d1, d2])


def cv_rf(X, y, n_splits=5, seed=0):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold
    from sklearn.metrics import r2_score
    kf = KFold(n_splits, shuffle=True, random_state=seed)
    oof = np.full(len(y), np.nan)
    for tr, te in kf.split(X):
        m = RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=0,
                                  min_samples_leaf=2)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    mask = ~np.isnan(oof)
    return (float(r2_score(y[mask], oof[mask])),
            float(np.mean(np.abs(y[mask] - oof[mask]))))


def main():
    df = load_primary()
    ee, dd, _ = F.monomer_features()
    pm = pair_medians(df)
    out = {}

    # baseline: full measurement-level (reference 0.386/0.422) is already known.
    # full pair-level median (all pairs, incl n=1)
    X_all = struct_feats(pm, ee, dd)
    for tgt in ["ln_r1", "ln_r2"]:
        r2, mae = cv_rf(X_all, pm[tgt].to_numpy(float))
        out[f"pairlevel_all_{tgt}"] = {"n_pairs": len(pm), "R2": round(r2, 3),
                                       "MAE": round(mae, 3)}

    for n_min in [2, 3, 4, 5]:
        sub = pm[pm.n_meas >= n_min].reset_index(drop=True)
        if len(sub) < 40:
            continue
        X = struct_feats(sub, ee, dd)
        entry = {"n_pairs": len(sub),
                 "n_meas_median": int(sub.n_meas.median()),
                 "total_var": round(float(sub.ln_r1.var()), 3),
                 "within_var_per_pair": round(
                     float((sub.ln_r1_std ** 2).mean()), 3)}
        # noise floor (median of n measurements -> within_var / n)
        entry["R2_ceiling"] = round(
            1 - float((sub.ln_r1_std ** 2 / sub.n_meas).mean()) /
            float(sub.ln_r1.var()), 3)
        for tgt in ["ln_r1", "ln_r2"]:
            r2, mae = cv_rf(X, sub[tgt].to_numpy(float))
            entry[f"{tgt}_R2"] = round(r2, 3)
            entry[f"{tgt}_MAE"] = round(mae, 3)
        out[f"gold_nmin{n_min}"] = entry
        print(f"n_min={n_min}: n_pairs={len(sub)} ceiling={entry['R2_ceiling']} "
              f"ln_r1 R2={entry['ln_r1_R2']} MAE={entry['ln_r1_MAE']} | "
              f"ln_r2 R2={entry['ln_r2_R2']} MAE={entry['ln_r2_MAE']}")

    print(json.dumps(out, indent=2))
    (P1 / "gold_subset.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
