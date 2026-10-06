"""P1 clean-subset experiment: primary-source-first (not replicate-count) filter.

Defines "clean" via provenance + methodology, not repeat count:
  primary  = source in {coxwell_full, coxwell_error, polycarp}   (exclude handbook/textbook/CKAs = secondary)
  modern   = fitting_method is a nonlinear / EVM / integrated / Tidwell-Mortimer /
             bootstrap / optimisation class (exclude Fineman-Ross / Kelen-Tudos /
             graphical linearisation)
  uncertain= r1_unc & r2_unc both present

Compares RF(struct+desc) R2/MAE on measurement-level across:
  full (5742), primary, primary+modern, primary+modern+uncertain,
  and the coxwell-only curated set.
Splits: random KFold(5), both_monomers_out, one_monomer_out.

Run: python3 src/p1_clean_subset.py
Output: projects/P1_rr_benchmark/clean_subset.json
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import STANDARDIZED, PROJECTS, RAW
import p1_features as F
from p1_benchmark import random_split, one_monomer_out, both_monomers_out

P1 = PROJECTS / "P1_rr_benchmark"
PRIMARY_SOURCES = {"coxwell_full", "coxwell_error", "polycarp"}


def coxwell_curated_ids():
    """source_row_id set of coxwell_full rows the Coxwell authors themselves
    marked as primary (excludes 'need to delete' + handbook/Young secondary)."""
    raw = pd.read_excel(RAW / "coxwell/Coxwell_RR full dataset final.xlsx")
    status = raw.iloc[:, 15].fillna("").astype(str).str.lower()
    src = raw["Source"].fillna("").astype(str).str.lower()
    drop = status.str.contains("need to delete")
    secondary = (status.str.contains("handbook") | status.str.contains("young")
                 | src.str.contains("polymer handbook")
                 | src.str.contains("brandrup")
                 | src.str.contains("young, l. j."))
    keep = ~(drop | secondary)
    return set(raw.index[keep].values), {"n_raw": len(raw),
                                         "n_drop_delete": int(drop.sum()),
                                         "n_secondary": int(secondary.sum())}

MODERN_KW = ["evm", "nlls", "nonlinear", "non-linear", "tidwell", "mortimer",
             "error-in-variables", "error in variables", "integrated",
             "copoint", "bootstrap", "marquardt", "simplex", "optimization",
             "optimisation", "least-squares", "least squares", "regression",
             "robust", "m-estimation"]
GRAPHICAL_KW = ["fineman", "finemann", "kelen", "tüdös", "tudos", "tudor",
                "graphical", "graphic", "intersection", "intersecting",
                "mayo-lewis plot", "plot", "slope"]


def load_primary():
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    return df[df.is_primary_measurement & df.ln_r1.notna() & df.ln_r2.notna()
              & (df.r1 > 0) & (df.r2 > 0)].reset_index(drop=True)


def struct_features(df):
    ee, dd, _ = F.monomer_features()
    e1 = ee.reindex(df.m1_inchikey).to_numpy(np.float32)
    e2 = ee.reindex(df.m2_inchikey).to_numpy(np.float32)
    d1 = dd.reindex(df.m1_inchikey).to_numpy(np.float32)
    d2 = dd.reindex(df.m2_inchikey).to_numpy(np.float32)
    return np.nan_to_num(np.hstack([e1, e2, d1, d2]), nan=0.0)


def add_filters(df):
    df = df.copy()
    m = df.fitting_method.fillna("").astype(str).str.lower()
    df["is_primary"] = df.source_dataset.isin(PRIMARY_SOURCES)
    df["is_modern"] = m.apply(lambda s: any(k in s for k in MODERN_KW)
                              and not any(g in s for g in GRAPHICAL_KW))
    df["has_unc"] = df.r1_unc.notna() & df.r2_unc.notna()
    return df


def eval_subset(df, X, target, splits, seed=0):
    from sklearn.ensemble import RandomForestRegressor
    y = df[target].to_numpy(float)
    oof = np.full(len(df), np.nan)
    for tr, te in splits:
        m = RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=seed,
                                  min_samples_leaf=2)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    mask = ~np.isnan(oof)
    ss = ((y[mask] - oof[mask]) ** 2).sum()
    tot = ((y[mask] - y[mask].mean()) ** 2).sum()
    return {"n_eval": int(mask.sum()),
            "MAE": float(np.mean(np.abs(y[mask] - oof[mask]))),
            "R2": float(1 - ss / tot)}


def main():
    df = add_filters(load_primary())
    X = struct_features(df)

    curated_ids, cur_meta = coxwell_curated_ids()
    is_curated_coxwell = (df.source_dataset == "coxwell_full") & \
        df.source_row_id.isin(curated_ids)

    subsets = {
        "full": np.ones(len(df), bool),
        "primary": df.is_primary.to_numpy(),
        "primary+modern": (df.is_primary & df.is_modern).to_numpy(),
        "primary+modern+unc": (df.is_primary & df.is_modern
                               & df.has_unc).to_numpy(),
        "coxwell_only": df.source_dataset.isin(
            {"coxwell_full", "coxwell_error"}).to_numpy(),
        "coxwell_curated": (is_curated_coxwell
                            | (df.source_dataset == "coxwell_error")).to_numpy(),
    }

    out = {}
    for tgt in ["ln_r1", "ln_r2"]:
        out[tgt] = {}
        for name, mask in subsets.items():
            idx = np.where(mask)[0]
            if len(idx) < 60:
                out[tgt][name] = {"n": int(len(idx)), "note": "too few"}
                continue
            sub = df.iloc[idx].reset_index(drop=True)
            Xs = X[idx]
            rec = {"n": int(len(idx)), "n_pairs": int(sub.canonical_pair_id.nunique())}
            for split_name, mk in [
                    ("random", lambda d: random_split(d)),
                    ("one_monomer_out", lambda d: one_monomer_out(d)),
                    ("both_monomers_out", lambda d: both_monomers_out(d))]:
                sp = mk(sub)
                if not sp:
                    rec[split_name] = {"note": "no folds"}
                    continue
                rec[split_name] = eval_subset(sub, Xs, tgt, sp)
            out[tgt][name] = rec
            print(f"[{tgt}] {name}: n={rec['n']} "
                  f"random R2={rec['random'].get('R2'):.3f} "
                  f"MAE={rec['random'].get('MAE'):.3f} "
                  f"| both-mono-out R2={rec['both_monomers_out'].get('R2', float('nan')):.3f}")

    (P1 / "clean_subset.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {P1/'clean_subset.json'}")


if __name__ == "__main__":
    main()
