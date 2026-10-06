"""P1 condition ablation: isolate temperature vs solvent contributions.

Run: python3 src/p1_ablation_conditions.py
Compares struct-only against +temperature / +solvent / +both on RF,
across random/pair/one-monomer-out/both-monomers-out, targets ln_r1/ln_r2.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import SPLITS, RESULTS
import p1_features as F
from p1_benchmark import (random_split, group_split, one_monomer_out,
                          both_monomers_out)

N_FOLDS = 5


def build_blocks(df):
    ee, dd, _ = F.monomer_features()
    struct = np.hstack([ee.reindex(df.m1_inchikey).to_numpy(np.float32),
                        ee.reindex(df.m2_inchikey).to_numpy(np.float32),
                        dd.reindex(df.m1_inchikey).to_numpy(np.float32),
                        dd.reindex(df.m2_inchikey).to_numpy(np.float32)])
    cond = F.condition_features(df)
    temp = cond[["temperature_K", "inv_temperature"]].to_numpy(np.float32)
    solv = cond[[c for c in cond.columns if c.startswith("solv_")]].to_numpy(np.float32)
    return struct, temp, solv


COMBOS = {"struct": ["s"], "+temp": ["s", "t"],
          "+solv": ["s", "v"], "+temp+solv": ["s", "t", "v"]}


def run():
    df = pd.read_parquet(SPLITS / "dataset_C_condition.parquet")
    df = df[df.ln_r1.notna() & df.ln_r2.notna()].reset_index(drop=True)
    struct, temp, solv = build_blocks(df)
    blocks = {"s": struct, "t": temp, "v": solv}

    splits = {"random": random_split(df),
              "pair": group_split(df, "canonical_pair_id"),
              "one_monomer_out": one_monomer_out(df),
              "both_monomers_out": both_monomers_out(df)}

    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score

    rows = []
    for target in ["ln_r1", "ln_r2"]:
        y = df[target].to_numpy(float)
        for sname, sp in splits.items():
            if not sp:
                continue
            for cname, keys in COMBOS.items():
                X = np.nan_to_num(np.hstack([blocks[k] for k in keys]),
                                  nan=0.0, posinf=0.0, neginf=0.0)
                oof = np.full(len(df), np.nan)
                for tr, te in sp:
                    m = RandomForestRegressor(n_estimators=300, n_jobs=-1,
                                              random_state=0, min_samples_leaf=2)
                    m.fit(X[tr], y[tr])
                    oof[te] = m.predict(X[te])
                mask = ~np.isnan(oof)
                rows.append({"target": target, "split": sname,
                             "features": cname, "n_eval": int(mask.sum()),
                             "MAE": float(np.mean(np.abs(y[mask] - oof[mask]))),
                             "R2": float(r2_score(y[mask], oof[mask]))})
                print(f"  {target}/{sname}/{cname}: MAE={rows[-1]['MAE']:.3f} "
                      f"R2={rows[-1]['R2']:.3f}")

    res = pd.DataFrame(rows)
    res.to_csv(RESULTS / "p1_condition_ablation.csv", index=False)

    print("\n=== MAE (struct | +temp | +solv | +temp+solv) ===")
    for target in ["ln_r1", "ln_r2"]:
        print(f"-- {target} --")
        p = res[res.target == target].pivot_table(
            index="split", columns="features", values="MAE")
        print(p.round(3).to_string())
    print(f"\nwrote {RESULTS/'p1_condition_ablation.csv'}")


if __name__ == "__main__":
    run()
