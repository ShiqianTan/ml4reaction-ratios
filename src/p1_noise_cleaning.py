"""P1 follow-up: (A) pair-level benchmark, (B) noise detection + cleaning.

(A) Pair-level: aggregate primary measurements to the canonical-pair median of
    ln r, then predict structure -> pair median with random + OOD splits.
    This is the scientifically meaningful target and the fair comparison vs the
    reference repo's pair-level R^2 = 0.62.

(B) Noise detection: three flags on *individual* primary measurements
      * within-pair MAD outlier: |ln r - pair median| > k * 1.4826 * MAD
      * quality_grade == 'D'
      * relative uncertainty r_unc/r > 0.5 (where available)
    Then evaluate whether removing flagged rows from TRAINING (fixed noisy test)
    and from BOTH train+test improves RF performance, at measurement level.

Run: python3 src/p1_noise_cleaning.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import STANDARDIZED, RESULTS, PROJECTS
import p1_features as F
from p1_benchmark import assign_monomer_folds, one_monomer_out, both_monomers_out

P1 = PROJECTS / "P1_rr_benchmark"
RANDOM_STATE = 0
N_FOLDS = 5


def load_primary():
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    d = df[df.is_primary_measurement & df.ln_r1.notna() & df.ln_r2.notna()]
    d = d[d.r1 > 0].reset_index(drop=True)
    return d


# --------------------------------------------------------------------------- #
# (A) pair-level
# --------------------------------------------------------------------------- #

def pair_level_df(d):
    g = d.groupby("canonical_pair_id")
    agg = g.agg(ln_r1=("ln_r1", "median"), ln_r2=("ln_r2", "median"),
                n_meas=("ln_r1", "size"),
                within_std_r1=("ln_r1", "std"),
                within_std_r2=("ln_r2", "std"),
                m1_inchikey=("m1_inchikey", "first"),
                m2_inchikey=("m2_inchikey", "first"),
                m1_name=("monomer_1_name", "first"),
                m2_name=("monomer_2_name", "first")).reset_index()
    return agg


def build_pair_features(pdf):
    ee, dd, _ = F.monomer_features()
    e1 = ee.reindex(pdf.m1_inchikey).to_numpy(np.float32)
    e2 = ee.reindex(pdf.m2_inchikey).to_numpy(np.float32)
    d1 = dd.reindex(pdf.m1_inchikey).to_numpy(np.float32)
    d2 = dd.reindex(pdf.m2_inchikey).to_numpy(np.float32)
    return np.nan_to_num(np.hstack([e1, e2, d1, d2]), nan=0.0)


def cv_eval(X, y, splits, make_model):
    from sklearn.base import clone
    oof = np.full(len(y), np.nan)
    for tr, te in splits:
        m = make_model()
        m.fit(X[tr], y[tr])
        oof[te] = m.predict(X[te])
    m = ~np.isnan(oof)
    resid = y[m] - oof[m]
    from sklearn.metrics import r2_score
    return dict(MAE=float(np.mean(np.abs(resid))),
                RMSE=float(np.sqrt(np.mean(resid ** 2))),
                R2=float(r2_score(y[m], oof[m])), n_eval=int(m.sum()))


def run_pair_level():
    d = load_primary()
    pdf = pair_level_df(d)
    X = build_pair_features(pdf)
    from sklearn.ensemble import RandomForestRegressor
    mk = lambda: RandomForestRegressor(n_estimators=400, n_jobs=-1,
                                       min_samples_leaf=2, random_state=0)
    out = {"n_pairs": len(pdf), "n_meas": len(d)}

    # random 5-fold
    from sklearn.model_selection import KFold
    kf = KFold(N_FOLDS, shuffle=True, random_state=0)
    rnd = list(kf.split(X))
    out["random"] = {**cv_eval(X, pdf.ln_r1.to_numpy(), rnd, mk),
                     "target": "ln_r1"}
    out["random_r2_target"] = {**cv_eval(X, pdf.ln_r2.to_numpy(), rnd, mk),
                               "target": "ln_r2"}

    # OOD: one-monomer-out + both-monomers-out
    oo = one_monomer_out(pdf)
    bo = both_monomers_out(pdf)
    out["one_monomer_out"] = {**cv_eval(X, pdf.ln_r1.to_numpy(), oo, mk),
                              "target": "ln_r1", "n_folds": len(oo)}
    out["both_monomers_out"] = {**cv_eval(X, pdf.ln_r1.to_numpy(), bo, mk),
                                "target": "ln_r1", "n_folds": len(bo)}
    out["one_monomer_out_r2"] = {**cv_eval(X, pdf.ln_r2.to_numpy(), oo, mk),
                                 "target": "ln_r2", "n_folds": len(oo)}
    out["both_monomers_out_r2"] = {**cv_eval(X, pdf.ln_r2.to_numpy(), bo, mk),
                                   "target": "ln_r2", "n_folds": len(bo)}

    # noise floor at pair level: within-pair std shrinks with n -> ceiling
    # total pair-level variance vs mean-within-pair variance / n_meas
    pdf.to_parquet(P1 / "pair_level_dataset.parquet", index=False)
    print("=== PAIR-LEVEL (RF struct+desc) ===")
    for k in ["random", "one_monomer_out", "both_monomers_out"]:
        print(f"  {k:18s} ln_r1 R2={out[k]['R2']:.3f} MAE={out[k]['MAE']:.3f}")
    return out


# --------------------------------------------------------------------------- #
# (B) noise detection + cleaning
# --------------------------------------------------------------------------- #

def add_noise_flags(d):
    """Return d with noise flag columns. Flags on individual measurements."""
    d = d.copy()
    d["flag_mad_outlier"] = False
    d["flag_low_quality"] = (d.quality_grade == "D").to_numpy()
    # relative uncertainty
    d["rel_unc_r1"] = d.r1_unc / d.r1
    d["flag_high_unc"] = d.rel_unc_r1 > 0.5
    # within-pair robust outlier (needs >=3 measurements)
    for t, col in [("ln_r1", "flag_mad_outlier"), ("ln_r2", "flag_mad_outlier2")]:
        d[col] = False
        g = d.groupby("canonical_pair_id")[t]
        cnt = g.transform("size")
        med = g.transform("median")
        mad = g.transform(lambda s: (s - s.median()).abs().median())
        z = (d[t] - med) / (1.4826 * mad.replace(0, np.nan))
        d[col] = d[col] | ((cnt >= 3) & (z.abs() > 3.0))
    d["flag_any_noise"] = (d.flag_mad_outlier | d.flag_low_quality |
                           d.flag_high_unc.fillna(False))
    return d


def run_cleaning():
    d = add_noise_flags(load_primary())
    ee, dd, _ = F.monomer_features()
    e1 = ee.reindex(d.m1_inchikey).to_numpy(np.float32)
    e2 = ee.reindex(d.m2_inchikey).to_numpy(np.float32)
    d1 = dd.reindex(d.m1_inchikey).to_numpy(np.float32)
    d2 = dd.reindex(d.m2_inchikey).to_numpy(np.float32)
    X = np.nan_to_num(np.hstack([e1, e2, d1, d2]), nan=0.0)
    y = d.ln_r1.to_numpy(float)

    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold

    counts = {
        "mad_outlier": int(d.flag_mad_outlier.sum()),
        "low_quality_D": int(d.flag_low_quality.sum()),
        "high_unc_r1": int(d.flag_high_unc.fillna(False).sum()),
        "any_noise": int(d.flag_any_noise.sum()),
        "n_total": len(d),
    }
    d.to_parquet(P1 / "noise_flagged_measurements.parquet", index=False)

    # fixed noisy TEST set (random 20%), clean only training
    kf = KFold(5, shuffle=True, random_state=0)
    # use first fold as fixed test for the "clean training" isolation experiment
    idx = np.arange(len(d))
    tr_full, te_full = next(iter(kf.split(idx)))

    schemes = {
        "full": np.ones(len(d), bool),
        "minus_MAD": ~d.flag_mad_outlier.to_numpy(),
        "minus_qualityD": ~d.flag_low_quality.to_numpy(),
        "minus_highunc": ~d.flag_high_unc.fillna(False).to_numpy(),
        "minus_any_noise": ~d.flag_any_noise.to_numpy(),
    }
    res = []
    for name, keep in schemes.items():
        # (a) train clean, test FIXED noisy (isolate training benefit)
        tr = tr_full[keep[tr_full]]
        te = te_full  # fixed noisy test
        m = RandomForestRegressor(n_estimators=400, n_jobs=-1,
                                  min_samples_leaf=2, random_state=0)
        m.fit(X[tr], y[tr])
        p = m.predict(X[te])
        resid = y[te] - p
        from sklearn.metrics import r2_score
        a = dict(scheme=name, n_train=int(len(tr)),
                 test_set="fixed_noisy",
                 MAE=float(np.mean(np.abs(resid))),
                 R2=float(r2_score(y[te], p)))
        # (b) train clean, test clean (deployment view)
        te2 = te_full[keep[te_full]]
        if len(te2) > 20:
            m2 = RandomForestRegressor(n_estimators=400, n_jobs=-1,
                                       min_samples_leaf=2, random_state=0)
            m2.fit(X[tr], y[tr])
            p2 = m2.predict(X[te2])
            r2 = y[te2] - p2
            a.update(n_test_clean=int(len(te2)), MAE_clean=float(np.mean(np.abs(r2))),
                     R2_clean=float(r2_score(y[te2], p2)))
        res.append(a)

    R = pd.DataFrame(res)
    R.to_csv(P1 / "noise_cleaning_eval.csv", index=False)
    print("\n=== NOISE CLEANING (measurement-level, RF struct+desc) ===")
    print(R.to_string(index=False))
    return {"counts": counts}


def main():
    out = {}
    out["pair_level"] = run_pair_level()
    out["cleaning"] = run_cleaning()
    (P1 / "p1_noise_pairlevel.json").write_text(json.dumps(out, indent=2,
                                                          default=str))
    print(f"\nwrote {P1/'p1_noise_pairlevel.json'}")


if __name__ == "__main__":
    main()
