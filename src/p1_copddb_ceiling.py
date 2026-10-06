"""P1 ceiling verification: can CopDDB first-principles Δ-barrier features break
the ECFP representation ceiling?

Hypothesis (P1_report conclusion #6): the binding constraint on ln-r prediction
is *feature representation* — ECFP cannot capture between-pair chemistry
(polarity Q-e, sterics, frontier orbitals). CopDDB supplies a first-principles
proxy: for r1 = k11/k12,  ln r1 ~ -ΔΔE‡/RT = -[E_barrier(R1,M2) - E_barrier(R1,M1)]/RT.

This script tests, on the CopDDB-covered subset, whether (a) the physical
relation holds, and (b) Δ-barrier features add predictive power over ECFP.

Run: python3 src/p1_copddb_ceiling.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
from scipy import stats

from common import STANDARDIZED, PROJECTS, R_GAS
import p1_features as F

P1 = PROJECTS / "P1_rr_benchmark"
HARTREE_KJ = 2625.5          # 1 Hartree = 2625.5 kJ/mol
T_REF = 298.15
SLOPE_THEORY = -HARTREE_KJ / (R_GAS * T_REF / 1000.0)   # ≈ -1059 per Hartree

DELTA_COLS = [c for c in ["copddb_d1_DE_barrier", "copddb_d2_DE_barrier",
                          "copddb_d1_DE_reaction", "copddb_d2_DE_reaction",
                          "copddb_d1_DE_SHgap", "copddb_d2_DE_SHgap",
                          "copddb_d1_DE_SLgap", "copddb_d2_DE_SLgap"]]


def load():
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    df = df[(df.is_primary_measurement) & (df.ln_r1.notna()) & (df.ln_r2.notna())
            & (df.r1 > 0) & (df.r2 > 0)].reset_index(drop=True)
    return df


def ecfp_block(df):
    ee, dd, _ = F.monomer_features()
    e1 = ee.reindex(df.m1_inchikey).to_numpy(np.float32)
    e2 = ee.reindex(df.m2_inchikey).to_numpy(np.float32)
    d1 = dd.reindex(df.m1_inchikey).to_numpy(np.float32)
    d2 = dd.reindex(df.m2_inchikey).to_numpy(np.float32)
    return np.nan_to_num(np.hstack([e1, e2, d1, d2]), nan=0.0)


def cv_rf(X, y, splits):
    from sklearn.ensemble import RandomForestRegressor
    oof = np.full(len(y), np.nan)
    for tr, te in splits:
        m = RandomForestRegressor(300, n_jobs=-1, random_state=0,
                                  min_samples_leaf=2).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    m = ~np.isnan(oof)
    r2 = 1 - ((y[m] - oof[m]) ** 2).sum() / ((y[m] - y[m].mean()) ** 2).sum()
    return r2, float(np.mean(np.abs(y[m] - oof[m])))


def cv_ridge(X, y, splits):
    from sklearn.linear_model import Ridge
    oof = np.full(len(y), np.nan)
    for tr, te in splits:
        m = Ridge(1.0).fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    m = ~np.isnan(oof)
    r2 = 1 - ((y[m] - oof[m]) ** 2).sum() / ((y[m] - y[m].mean()) ** 2).sum()
    return r2, float(np.mean(np.abs(y[m] - oof[m])))


def main():
    df = load()
    phys = F.copddb_rr_descriptors(df)
    both = phys[["copddb_d1_DE_barrier", "copddb_d2_DE_barrier"]].notna().all(1)
    sub = df[both].reset_index(drop=True)
    D = phys.loc[both, DELTA_COLS].to_numpy(np.float32)
    X_ecfp = ecfp_block(sub)

    from sklearn.model_selection import KFold, GroupKFold
    splits_rand = list(KFold(5, shuffle=True, random_state=0).split(sub))
    splits_pair = list(GroupKFold(5).split(
        sub, groups=sub.canonical_pair_id.fillna("NA").astype(str).values))

    out = {"n_subset": int(len(sub)), "n_pairs": int(sub.canonical_pair_id.nunique()),
           "n_monomers": int(pd.concat([sub.m1_inchikey, sub.m2_inchikey]).nunique())}

    for target, dcol in [("ln_r1", "copddb_d1_DE_barrier"),
                         ("ln_r2", "copddb_d2_DE_barrier")]:
        y = sub[target].to_numpy(float)
        d = phys.loc[both, dcol].to_numpy(float)
        # physical relation: ln r vs Δ-barrier
        pear, p = stats.pearsonr(d, y)
        slope, intercept, rv, pv, se = stats.linregress(d, y)
        out[target + "_phys"] = {
            "pearson": round(float(pear), 3), "pearson_p": round(float(p), 4),
            "slope": round(float(slope), 1), "slope_theory": round(SLOPE_THEORY, 1),
            "r2_phys": round(float(rv ** 2), 3)}

        # models on subset (random + pair split)
        D1 = D[:, 0:2] if target == "ln_r1" else D[:, 0:2]
        out[target + "_ecfp_rf_rand"] = cv_rf(X_ecfp, y, splits_rand)
        out[target + "_delta_ridge_rand"] = cv_ridge(np.nan_to_num(D), y, splits_rand)
        out[target + "_ecfp_delta_rf_rand"] = cv_rf(np.hstack([X_ecfp, np.nan_to_num(D)]), y, splits_rand)
        out[target + "_ecfp_rf_pair"] = cv_rf(X_ecfp, y, splits_pair)
        out[target + "_delta_ridge_pair"] = cv_ridge(np.nan_to_num(D), y, splits_pair)
        out[target + "_ecfp_delta_rf_pair"] = cv_rf(np.hstack([X_ecfp, np.nan_to_num(D)]), y, splits_pair)

    print(json.dumps(out, indent=2))
    (P1 / "copddb_ceiling.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {P1/'copddb_ceiling.json'}")


if __name__ == "__main__":
    main()
