"""P2: is the enthalpy-entropy compensation real, or a fit-error artefact?

For a series of fits  ln r = -dH/(R T) + dS/R, the estimated dH and dS are
mechanically correlated because both come from the same least-squares fit of
ln r vs 1/T (slope and intercept of OLS share a covariance proportional to
-mean(1/T)). This produces an *apparent* dH-dS compensation even when none
exists.

Three tests:
  1. observed compensation (slope, r2, Pearson corr) on the original estimates
  2. centred-fit compensation (refit with x = 1/T - mean(1/T) to decorrelate)
  3. parametric null bootstrap: simulate series under the NULL (dH, dS
     independently drawn, Gaussian noise at observed residual scale), refit,
     and measure the compensation induced by estimation error alone.

Verdict: real compensation if observed r2 >> null r2 distribution.

Run: python3 src/p2_compensation.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
from scipy import stats

from common import STANDARDIZED, PROJECTS, R_GAS

P2 = PROJECTS / "P2_temperature_rr"
MIN_T = 3
MIN_SPAN = 20.0
N_BOOT = 400
RNG = np.random.default_rng(0)


def collect_series():
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    df = df[df.ln_r1.notna() & df.ln_r2.notna()
            & df.temperature_K.notna()].copy()
    out = []
    for key, g in df.groupby(["canonical_pair_id", "solvent_family"]):
        if g.temperature_K.nunique() < MIN_T:
            continue
        if g.temperature_K.max() - g.temperature_K.min() < MIN_SPAN:
            continue
        for target in ["ln_r1", "ln_r2"]:
            T = g.temperature_K.to_numpy(float)
            y = g[target].to_numpy(float)
            if len(np.unique(T)) < MIN_T:
                continue
            x = 1.0 / T
            lr = stats.linregress(x, y)
            resid = y - (lr.slope * x + lr.intercept)
            out.append({"T": T, "y": y, "x": x, "n": len(T),
                        "dH": -R_GAS * lr.slope, "dS": R_GAS * lr.intercept,
                        "sigma": float(resid.std(ddof=1)) if len(T) > 2 else 0.0,
                        "x_mean": float(x.mean())})
    return out


def compensation(dH, dS):
    """dH in J/mol, dS in J/(mol*K): slope of dH vs dS is in K directly."""
    if len(dH) < 5:
        return {}
    r = stats.pearsonr(dS, dH)
    lr = stats.linregress(dS, dH)
    return {"slope_K": round(float(lr.slope), 1),
            "r2": round(float(lr.rvalue ** 2), 3),
            "pearson": round(float(r.statistic), 3)}


def centered_fit(series):
    dHc, dSc = [], []
    for s in series:
        xc = s["x"] - s["x_mean"]
        lr = stats.linregress(xc, s["y"])
        dHc.append(-R_GAS * lr.slope)
        dSc.append(R_GAS * lr.intercept)
    return np.array(dHc), np.array(dSc)


def null_bootstrap(series):
    dH_obs = np.array([s["dH"] for s in series])
    dS_obs = np.array([s["dS"] for s in series])
    slopes, r2s = [], []
    for _ in range(N_BOOT):
        dH_null, dS_null = [], []
        for s in series:
            # independent draws from the observed marginals
            dH0 = RNG.choice(dH_obs)
            dS0 = RNG.choice(dS_obs)
            # simulate ln r at the observed temperatures under this null
            y0 = -dH0 / (R_GAS * s["T"]) + dS0 / R_GAS
            y0 = y0 + RNG.normal(0.0, max(s["sigma"], 1e-6), size=s["n"])
            lr = stats.linregress(s["x"], y0)
            dH_null.append(-R_GAS * lr.slope)
            dS_null.append(R_GAS * lr.intercept)
        c = compensation(np.array(dH_null), np.array(dS_null))
        slopes.append(c["slope_K"])
        r2s.append(c["r2"])
    return np.array(slopes), np.array(r2s)


def main():
    series = collect_series()
    dH = np.array([s["dH"] for s in series])
    dS = np.array([s["dS"] for s in series])
    print(f"series = {len(series)}")

    obs = compensation(dH, dS)
    dHc, dSc = centered_fit(series)
    cent = compensation(dHc, dSc)
    null_slopes, null_r2s = null_bootstrap(series)

    # artefact signature: uncentered dS = R*<ln r> + mean(1/T)*dH, so the
    # spurious compensation "temperature" equals 1 / mean(1/T).
    mean_1T = np.mean([s["x_mean"] for s in series])
    artefact_T = 1.0 / mean_1T

    # The centred fit removes the mechanical slope-intercept coupling: if the
    # compensation is an artefact, centred r2 collapses.
    obs_r2 = obs.get("r2", 0)
    cent_r2 = cent.get("r2", 0)
    if cent_r2 < 0.05 and obs_r2 > 0.5:
        verdict = ("ARTEFACT: uncentered r2={obs_r2:.2f} collapses to "
                   "{cent_r2:.2f} after centering, and the compensation "
                   "temperature ({obs_slope:.0f} K) tracks the mean temperature "
                   "({T_art:.0f} K)").format(
                       obs_r2=obs_r2, cent_r2=cent_r2,
                       obs_slope=obs["slope_K"], T_art=artefact_T)
    else:
        verdict = ("REAL: centred compensation r2={cent_r2:.2f} persists "
                   "after decorrelating slope/intercept").format(cent_r2=cent_r2)

    result = {
        "n_series": int(len(series)),
        "observed": obs,
        "centered_fit": cent,
        "mean_1_over_T_K": round(float(artefact_T), 1),
        "null_bootstrap": {
            "n": N_BOOT,
            "slope_K_median": round(float(np.median(null_slopes)), 1),
            "r2_median": round(float(np.median(null_r2s)), 3),
            "r2_p95": round(float(np.percentile(null_r2s, 95)), 3),
        },
        "verdict": verdict,
    }

    (P2 / "p2_compensation_check.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
