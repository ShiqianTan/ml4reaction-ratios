"""P1 ReactivityRatioBench: model x feature x split experiment matrix.

Run: python3 src/p1_benchmark.py [--quick]
Outputs: results/p1_benchmark_results.csv, results/p1_split_sizes.csv
"""
from __future__ import annotations

import sys
import time
import numpy as np
import pandas as pd

from common import STANDARDIZED, SPLITS, RESULTS, PROJECTS, R_GAS
import p1_features as F

P1 = PROJECTS / "P1_rr_benchmark"
QUICK = "--quick" in sys.argv
N_FOLDS = 5


# --------------------------------------------------------------------------- #
# splits
# --------------------------------------------------------------------------- #

def random_split(df, seed=0):
    from sklearn.model_selection import KFold
    kf = KFold(N_FOLDS, shuffle=True, random_state=seed)
    return list(kf.split(df))


def group_split(df, group_col, seed=0):
    from sklearn.model_selection import GroupKFold
    gkf = GroupKFold(N_FOLDS)
    groups = df[group_col].fillna("NA").values
    return list(gkf.split(df, groups=groups))


def assign_monomer_folds(monomers, edges, n_folds, seed=0):
    """Greedy balanced graph partition of monomers into folds."""
    rng = np.random.default_rng(seed)
    deg = {m: 0 for m in monomers}
    for a, b in edges:
        if a in deg:
            deg[a] += 1
        if b in deg:
            deg[b] += 1
    order = sorted(monomers, key=lambda m: (-deg[m], rng.random()))
    fold_size = np.zeros(n_folds)
    assign = {}
    for m in order:
        k = int(np.argmin(fold_size))
        assign[m] = k
        fold_size[k] += 1
    return assign


def both_monomers_out(df, seed=0):
    monomers = pd.unique(pd.concat([df.m1_inchikey, df.m2_inchikey]).dropna())
    edges = list(zip(df.m1_inchikey, df.m2_inchikey))
    assign = assign_monomer_folds(monomers, edges, N_FOLDS, seed)
    splits = []
    for k in range(N_FOLDS):
        test = df.index[[(assign.get(a, -1) == k) and (assign.get(b, -1) == k)
                         for a, b in zip(df.m1_inchikey, df.m2_inchikey)]]
        train = df.index[[(assign.get(a, -1) != k) and (assign.get(b, -1) != k)
                         for a, b in zip(df.m1_inchikey, df.m2_inchikey)]]
        if len(test) >= 10 and len(train) >= 50:
            splits.append((np.array(train), np.array(test)))
    return splits


def one_monomer_out(df, seed=0):
    monomers = pd.unique(pd.concat([df.m1_inchikey, df.m2_inchikey]).dropna())
    edges = list(zip(df.m1_inchikey, df.m2_inchikey))
    assign = assign_monomer_folds(monomers, edges, N_FOLDS, seed)
    splits = []
    for k in range(N_FOLDS):
        test = df.index[[(assign.get(a, -1) == k) or (assign.get(b, -1) == k)
                         for a, b in zip(df.m1_inchikey, df.m2_inchikey)]]
        train = df.index[[(assign.get(a, -1) != k) and (assign.get(b, -1) != k)
                          for a, b in zip(df.m1_inchikey, df.m2_inchikey)]]
        if len(test) >= 10 and len(train) >= 50:
            splits.append((np.array(train), np.array(test)))
    return splits


# --------------------------------------------------------------------------- #
# models
# --------------------------------------------------------------------------- #

def make_models():
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    models = {
        "ridge": Ridge(alpha=1.0),
        "rf": RandomForestRegressor(n_estimators=300, n_jobs=-1,
                                    random_state=0, min_samples_leaf=2),
        "xgb": None,
    }
    try:
        import xgboost as xgb
        models = {"ridge": models["ridge"], "rf": models["rf"],
                  "xgb": xgb.XGBRegressor(
                      n_estimators=300, learning_rate=0.06, max_depth=6,
                      subsample=0.8, colsample_bytree=0.3, n_jobs=-1,
                      tree_method="hist", random_state=0)}
    except Exception as e:  # pragma: no cover
        print("xgboost unavailable:", e)
        models.pop("xgb")
    if QUICK:
        models.pop("xgb", None)
        models["rf"].set_params(n_estimators=80)
    else:
        from sklearn.neural_network import MLPRegressor
        models["mlp"] = MLPRegressor(
            hidden_layer_sizes=(256,), alpha=1e-3, max_iter=200,
            early_stopping=True, n_iter_no_change=20, random_state=0)
    return models


def metrics(y, p):
    from sklearn.metrics import r2_score
    resid = y - p
    sigma = float(np.std(resid))
    nll = 0.5 * np.log(2 * np.pi * sigma ** 2) + float(np.mean(resid ** 2)) / (2 * sigma ** 2)
    return {"MAE": float(np.mean(np.abs(resid))),
            "RMSE": float(np.sqrt(np.mean(resid ** 2))),
            "R2": float(r2_score(y, p)) if len(y) > 2 else np.nan,
            "NLL": float(nll) if sigma > 1e-9 else np.nan,
            "coverage95": float(np.mean(np.abs(resid) < 1.96 * sigma))
            if sigma > 1e-9 else np.nan}


def mayo_lewis_curve_error(r1t, r2t, r1p, r2p):
    f1 = np.linspace(0.02, 0.98, 49)[None, :]
    f2 = 1 - f1

    def F1(r1, r2):
        r1 = np.asarray(r1, dtype=float)[:, None]
        r2 = np.asarray(r2, dtype=float)[:, None]
        return (r1 * f1 ** 2 + f1 * f2) / (r1 * f1 ** 2 + 2 * f1 * f2
                                           + r2 * f2 ** 2)
    return float(np.mean(np.abs(F1(r1t, r2t) - F1(r1p, r2p))))


# --------------------------------------------------------------------------- #
# experiment driver
# --------------------------------------------------------------------------- #

def build_matrix(df):
    ee, dd, smi = F.monomer_features()
    ecols = list(ee.columns)
    dcols = list(dd.columns)
    ee1 = ee.reindex(df.m1_inchikey).to_numpy(np.float32)
    ee2 = ee.reindex(df.m2_inchikey).to_numpy(np.float32)
    dd1 = dd.reindex(df.m1_inchikey).to_numpy(np.float32)
    dd2 = dd.reindex(df.m2_inchikey).to_numpy(np.float32)
    phys = F.copddb_rr_descriptors(df)
    avail = phys.pop("copddb_available").to_numpy(np.float32)[:, None]
    blocks = {
        "struct": np.hstack([ee1, ee2]),
        "desc": np.hstack([dd1, dd2]),
        "cond": F.condition_features(df).to_numpy(np.float32),
        "phys": phys.to_numpy(np.float32),
        "phys_avail": avail,
    }
    return blocks


FEATURE_SETS = {
    "struct": ["struct"],
    "struct+desc": ["struct", "desc"],
    "struct+cond": ["struct", "desc", "cond"],
    "struct+phys": ["struct", "desc", "phys", "phys_avail"],
    "struct+cond+phys": ["struct", "desc", "cond", "phys", "phys_avail"],
}
# focused matrix for the full (GBDT-inclusive) run
FULL_FEATURE_SETS = {k: FEATURE_SETS[k]
                     for k in ["struct", "struct+cond", "struct+cond+phys"]}
if not QUICK:
    FEATURE_SETS = FULL_FEATURE_SETS


def run(df, target, split_name, splits, blocks, models):
    rows = []
    oofs = {}
    y_all = df[target].to_numpy(float)
    for fs_name, keys in FEATURE_SETS.items():
        X = np.hstack([blocks[k] for k in keys])
        for mname, model in models.items():
            if mname == "mlp" and fs_name != "struct":
                continue
            oof = np.full(len(df), np.nan)
            for tr, te in splits:
                from sklearn.base import clone
                m = clone(model)
                Xtr = np.nan_to_num(X[tr], nan=0.0, posinf=0.0, neginf=0.0)
                Xte = np.nan_to_num(X[te], nan=0.0, posinf=0.0, neginf=0.0)
                m.fit(Xtr, y_all[tr])
                oof[te] = m.predict(Xte)
            mask = ~np.isnan(oof)
            met = metrics(y_all[mask], oof[mask])
            oofs[(fs_name, mname)] = oof
            rows.append({"target": target, "split": split_name,
                         "features": fs_name, "model": mname,
                         "n_eval": int(mask.sum()), **met})
    return rows, oofs


def main():
    t0 = time.time()
    df = pd.read_parquet(SPLITS / "dataset_C_condition.parquet")
    df = df.reset_index(drop=True)
    # targets: use log/ln r (drop non-positive)
    df = df[df.ln_r1.notna() & df.ln_r2.notna()].reset_index(drop=True)
    print(f"benchmark rows: {len(df):,}")

    blocks = build_matrix(df)
    models = make_models()

    split_defs = {
        "random": random_split(df),
        "pair": group_split(df, "canonical_pair_id"),
        "one_monomer_out": one_monomer_out(df),
        "both_monomers_out": both_monomers_out(df),
    }
    sizes = {k: [(len(a), len(b)) for a, b in v] for k, v in split_defs.items()}
    print("split folds (train,test):")
    for k, v in sizes.items():
        print(" ", k, v)
    pd.DataFrame([{"split": k, "fold": i, "n_train": a, "n_test": b}
                  for k, v in sizes.items() for i, (a, b) in enumerate(v)]
                 ).to_csv(RESULTS / "p1_split_sizes.csv", index=False)

    all_rows = []
    oof_store = {}  # (split, features, model) -> {target: oof_array}
    for target in ["ln_r1", "ln_r2"]:
        for sname, splits in split_defs.items():
            if not splits:
                continue
            rows, oofs = run(df, target, sname, splits, blocks, models)
            all_rows += rows
            for (fs, mn), oof in oofs.items():
                oof_store.setdefault((sname, fs, mn), {})[target] = oof
            print(f"  done {target} / {sname}  ({time.time()-t0:.0f}s)")

    res = pd.DataFrame(all_rows)

    # ---- Mayo-Lewis curve error + joint OOF table (needs both directions) ----
    joint_rows = []
    ecurve = {}
    r1_true = df["r1"].to_numpy(float)
    r2_true = df["r2"].to_numpy(float)
    for (sname, fs, mn), d in oof_store.items():
        if "ln_r1" not in d or "ln_r2" not in d:
            continue
        o1, o2 = d["ln_r1"], d["ln_r2"]
        m = ~np.isnan(o1) & ~np.isnan(o2)
        r1p = np.exp(o1[m])
        r2p = np.exp(o2[m])
        e = mayo_lewis_curve_error(r1_true[m], r2_true[m], r1p, r2p)
        ecurve[(sname, fs, mn)] = e
        sub = df[m].copy()
        sub["split"], sub["features"], sub["model"] = sname, fs, mn
        sub["ln_r1_oof"], sub["ln_r2_oof"] = o1[m], o2[m]
        sub["r1_oof"], sub["r2_oof"] = r1p, r2p
        joint_rows.append(sub[[
            "canonical_pair_id", "m1_inchikey", "m2_inchikey", "r1", "r2",
            "ln_r1", "ln_r2", "temperature_K", "solvent_family",
            "split", "features", "model", "ln_r1_oof", "ln_r2_oof",
            "r1_oof", "r2_oof"]])
    res["E_curve"] = [ecurve.get((s, f, m), np.nan)
                      for s, f, m in zip(res.split, res.features, res.model)]
    if joint_rows:
        pd.concat(joint_rows).to_parquet(RESULTS / "p1_benchmark_oof.parquet",
                                         index=False)

    res.to_csv(RESULTS / "p1_benchmark_results.csv", index=False)
    print("\n=== summary (MAE by split, best feature set per model/target) ===")
    piv = res.pivot_table(index=["target", "model", "split"],
                          columns="features", values="MAE")
    print(piv.round(3).to_string())
    print(f"\nwrote {RESULTS/'p1_benchmark_results.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
