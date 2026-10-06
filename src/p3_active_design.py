"""P3: ML prior + Bayesian / active experimental design — strict offline benchmark.

Upgrades over the first version:
  * fixed shared cold-start per seed (fair strategy comparison)
  * multiple seeds -> mean +/- std of final MAE
  * calibration metrics (NLL, coverage95) on the held-out pairs
  * greedy-oracle upper bound (labels the pool points with largest true error)
  * normalised improvement (MAE gain / oracle-achievable gain)

Run: python3 src/p3_active_design.py
Outputs under projects/P3_active_rr_design/
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import SPLITS, PROJECTS
import p1_features as F

P3 = PROJECTS / "P3_active_rr_design"
N_INIT = 120
BATCH = 30
N_ROUNDS = 12          # N_INIT + BATCH*N_ROUNDS = 480 experiments
HOLDOUT_PAIR_FRAC = 0.25
TARGET = "ln_r1"
SEEDS = [0, 1, 2]

STRATEGIES = ["random", "uncertainty", "diversity", "unc+div", "physics-prior",
              "greedy-oracle"]


def get_data():
    df = pd.read_parquet(SPLITS / "dataset_C_condition.parquet")
    df = df[df.ln_r1.notna() & df.ln_r2.notna()].reset_index(drop=True)
    return df


def features(df):
    ee, dd, _ = F.monomer_features()
    e1 = ee.reindex(df.m1_inchikey).to_numpy(np.float32)
    e2 = ee.reindex(df.m2_inchikey).to_numpy(np.float32)
    d1 = dd.reindex(df.m1_inchikey).to_numpy(np.float32)
    d2 = dd.reindex(df.m2_inchikey).to_numpy(np.float32)
    cond = F.condition_features(df).to_numpy(np.float32)
    phys = F.copddb_rr_descriptors(df)
    has_phys = phys["copddb_available"].to_numpy(bool)
    X = np.nan_to_num(np.hstack([e1, e2, d1, d2, cond]), nan=0.0)
    return X, has_phys


def ensemble():
    from sklearn.ensemble import RandomForestRegressor
    return [RandomForestRegressor(n_estimators=100, min_samples_leaf=2,
                                  n_jobs=-1, random_state=s)
            for s in range(5)]


def fit_predict(models, Xl, yl, Xq):
    for m in models:
        m.fit(Xl, yl)
    P = np.array([m.predict(Xq) for m in models])
    return P.mean(0), P.std(0)


def tanimoto_min(F_q, F_l):
    A = F_q.astype(bool)
    B = F_l.astype(bool)
    inter = A @ B.T
    union = A.sum(1)[:, None] + B.sum(1)[None, :] - inter
    sim = inter / np.maximum(union, 1)
    return sim.max(1)


def calibration(y, pred):
    resid = y - pred
    sigma = float(resid.std() + 1e-9)
    nll = 0.5 * np.log(2 * np.pi * sigma ** 2) + float(np.mean(resid ** 2)) / (2 * sigma ** 2)
    cov = float(np.mean(np.abs(resid) < 1.96 * sigma))
    return nll, cov


def acquire(strategy, cand, cand_pos, ecfp_cand, ecfp_labs, std_q, has_phys,
            resid_q, rng, BATCH):
    n = len(cand)
    if strategy == "random":
        sel = rng.choice(n, BATCH, replace=False)
    elif strategy == "uncertainty":
        sel = np.argsort(-std_q)[:BATCH]
    elif strategy == "diversity":
        sim = tanimoto_min(ecfp_cand, ecfp_labs)
        sel = np.argsort(sim)[:BATCH]
    elif strategy == "unc+div":
        u = (std_q - std_q.min()) / (np.ptp(std_q) + 1e-9)
        sim = tanimoto_min(ecfp_cand, ecfp_labs)
        d = 1 - (sim - sim.min()) / (np.ptp(sim) + 1e-9)
        sel = np.argsort(-(0.5 * u + 0.5 * d))[:BATCH]
    elif strategy == "physics-prior":
        pr = has_phys[cand_pos].astype(float)
        u = (std_q - std_q.min()) / (np.ptp(std_q) + 1e-9)
        sel = np.argsort(-(0.5 * pr + 0.5 * u))[:BATCH]
    elif strategy == "greedy-oracle":
        sel = np.argsort(-np.abs(resid_q))[:BATCH]
    return cand[sel]


def run_seed(seed, df, X, has_phys, y, test_idx, pool_idx, ecfp):
    rng = np.random.default_rng(seed)
    cold = list(rng.choice(pool_idx, N_INIT, replace=False))
    rows = []
    for strategy in STRATEGIES:
        labs = list(cold)
        cand = np.setdiff1d(pool_idx, np.array(labs))
        for rnd in range(N_ROUNDS + 1):
            models = ensemble()
            mean_t, _ = fit_predict(models, X[labs], y[labs], X[test_idx])
            mean_q, std_q = fit_predict(models, X[labs], y[labs], X[cand])
            mae = float(np.mean(np.abs(y[test_idx] - mean_t)))
            nll, cov = calibration(y[test_idx], mean_t)
            rows.append({"seed": seed, "strategy": strategy, "round": rnd,
                         "n_exp": len(labs), "MAE": mae, "NLL": nll,
                         "coverage95": cov})
            if rnd == N_ROUNDS:
                break
            resid_q = y[cand] - mean_q
            chosen = acquire(strategy, cand, cand, ecfp[cand], ecfp[labs],
                             std_q, has_phys, resid_q, rng, BATCH)
            labs.extend(chosen.tolist())
            cand = np.setdiff1d(cand, chosen)
    return pd.DataFrame(rows)


def main():
    df = get_data()
    X, has_phys = features(df)
    y = df[TARGET].to_numpy(float)
    pairs = df.canonical_pair_id.fillna("NA").astype(str).values
    uniq = np.array(sorted(set(pairs)))
    ecfp = X[:, :4096]

    all_rows = []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        hold = set(rng.choice(uniq, int(len(uniq) * HOLDOUT_PAIR_FRAC),
                              replace=False))
        test_mask = np.array([p in hold for p in pairs])
        test_idx = np.where(test_mask)[0]
        pool_idx = np.where(~test_mask)[0]
        print(f"seed={seed} pairs={len(uniq)} test={len(test_idx)} "
              f"pool={len(pool_idx)}")
        all_rows.append(run_seed(seed, df, X, has_phys, y, test_idx, pool_idx,
                                 ecfp))

    res = pd.concat(all_rows, ignore_index=True)
    res.to_csv(P3 / "active_learning_curve.csv", index=False)

    fin = res[res.n_exp == res.n_exp.max()]
    ini = res[res.n_exp == N_INIT]
    summary = {"n_init": N_INIT, "batch": BATCH, "n_rounds": N_ROUNDS,
               "n_seeds": len(SEEDS), "n_test_rows": int(test_idx.size)}
    for s in STRATEGIES:
        f = fin[fin.strategy == s]
        i = ini[ini.strategy == s]
        oracle = fin[fin.strategy == "greedy-oracle"].MAE.mean()
        summary[s] = {
            "MAE_init_mean": round(float(i.MAE.mean()), 3),
            "MAE_final_mean": round(float(f.MAE.mean()), 3),
            "MAE_final_std": round(float(f.MAE.std()), 3),
            "improvement": round(float(i.MAE.mean() - f.MAE.mean()), 3),
            "NLL_final": round(float(f.NLL.mean()), 3),
            "coverage95_final": round(float(f.coverage95.mean()), 3),
        }
        if s != "greedy-oracle":
            gain_oracle = max(i.MAE.mean() - oracle, 1e-9)
            summary[s]["normalised_improvement"] = round(
                (i.MAE.mean() - f.MAE.mean()) / gain_oracle, 3)

    print(json.dumps(summary, indent=2))
    (P3 / "p3_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {P3/'active_learning_curve.csv'} and p3_summary.json")


if __name__ == "__main__":
    main()
