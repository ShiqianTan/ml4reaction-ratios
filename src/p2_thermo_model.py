"""P2b: can the temperature response (dH, dS) be predicted from structure?

Leave-pair-out CV predicting the fitted dH (kJ/mol) and dS (J/mol/K)
from ECFP4 of the radical/monomer in the reaction direction.

Run: python3 src/p2_thermo_model.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import STANDARDIZED, PROJECTS
import p1_features as F

P2 = PROJECTS / "P2_temperature_rr"


def build():
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    series = pd.read_csv(P2 / "temperature_series_fits.csv")
    strict = series[series.grouping == "strict"].copy()
    # attach inchikeys
    key = df.drop_duplicates("canonical_pair_id").set_index("canonical_pair_id")[
        ["m1_inchikey", "m2_inchikey"]]
    strict = strict.join(key, on="pair_id")
    ee, dd, smi = F.monomer_features()
    rows, meta = [], []
    for _, r in strict.iterrows():
        if r.target == "ln_r1":
            rad, com = r.m1_inchikey, r.m2_inchikey
        else:
            rad, com = r.m2_inchikey, r.m1_inchikey
        if rad not in ee.index or com not in ee.index:
            continue
        feat = np.concatenate([ee.loc[rad].to_numpy(np.float32),
                               ee.loc[com].to_numpy(np.float32),
                               dd.loc[rad].to_numpy(np.float32),
                               dd.loc[com].to_numpy(np.float32)])
        rows.append(feat)
        meta.append({"pair_id": r.pair_id, "target": r.target,
                     "dH": r.dH_kJ_mol, "dS": r.dS_J_molK,
                     "mean_T": r.mean_T, "p_value": r.p_value})
    return np.array(rows, np.float32), pd.DataFrame(meta)


def cv_predict(X, y, groups, model):
    from sklearn.base import clone
    from sklearn.model_selection import LeaveOneGroupOut
    oof = np.full(len(y), np.nan)
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        if len(tr) < 10:
            continue
        m = clone(model)
        m.fit(np.nan_to_num(X[tr]), y[tr])
        oof[te] = m.predict(np.nan_to_num(X[te]))
    return oof


def main():
    from sklearn.linear_model import Ridge
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    X, meta = build()
    print(f"thermo-model samples: {len(X)}  pairs: {meta.pair_id.nunique()}")
    out = {}
    for target in ["dH", "dS"]:
        y = meta[target].to_numpy(float)
        for name, model in [("ridge", Ridge(alpha=10.0)),
                            ("rf", RandomForestRegressor(
                                n_estimators=300, min_samples_leaf=3,
                                n_jobs=-1, random_state=0))]:
            oof = cv_predict(X, y, meta.pair_id.values, model)
            mask = ~np.isnan(oof)
            out[f"{target}_{name}"] = {
                "n": int(mask.sum()),
                "MAE": round(float(np.mean(np.abs(y[mask] - oof[mask]))), 3),
                "R2": round(float(r2_score(y[mask], oof[mask])), 3),
                "MAE_vs_mean": round(float(np.mean(np.abs(
                    y[mask] - y[mask].mean()))), 3),
            }
            print(target, name, out[f"{target}_{name}"])
    # significant-only subset (genuine temperature dependence)
    sig = meta.p_value < 0.05
    if sig.sum() >= 30:
        for target in ["dH", "dS"]:
            y = meta[target].to_numpy(float)[sig.values]
            Xs = X[sig.values]
            g = meta.pair_id.values[sig.values]
            oof = cv_predict(Xs, y, g, Ridge(alpha=10.0))
            m = ~np.isnan(oof)
            out[f"{target}_ridge_significant_only"] = {
                "n": int(m.sum()),
                "R2": round(float(r2_score(y[m], oof[m])), 3),
                "MAE_vs_mean": round(float(np.mean(np.abs(
                    y[m] - y[m].mean()))), 3),
            }
    # Is the *sign* of dH predictable? (relevant for whether T raises/lowers r)
    y = np.sign(meta.dH.to_numpy(float))
    oof = cv_predict(X, meta.dH.to_numpy(float), meta.pair_id.values,
                     Ridge(alpha=10.0))
    acc = float(np.mean(np.sign(oof[~np.isnan(oof)])
                        == y[~np.isnan(oof)]))
    out["dH_sign_accuracy_ridge"] = round(acc, 3)
    print("dH sign accuracy (ridge):", round(acc, 3))
    (P2 / "p2_thermo_model.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
