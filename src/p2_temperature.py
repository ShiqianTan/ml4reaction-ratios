"""P2: temperature-dependent reactivity-ratio thermodynamics.

  ln r = -dH/(R T) + dS/R      per (pair, solvent-family) series
  dH = -R * slope ;  dS = R * intercept

Run: python3 src/p2_temperature.py
Outputs under projects/P2_temperature_rr/
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
from scipy import stats

from common import STANDARDIZED, RESULTS, PROJECTS, R_GAS

P2 = PROJECTS / "P2_temperature_rr"
MIN_T = 3
MIN_SPAN = 20.0


def load():
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    df = df[df.ln_r1.notna() & df.ln_r2.notna()
            & df.temperature_K.notna()].copy()
    return df


def fit_series(g, target):
    """Return dict of fit stats for ln r vs 1/T."""
    T = g.temperature_K.to_numpy(float)
    y = g[target].to_numpy(float)
    if len(np.unique(T)) < MIN_T:
        return None
    x = 1.0 / T
    lr = stats.linregress(x, y)
    dH = -R_GAS * lr.slope            # J/mol
    dS = R_GAS * lr.intercept         # J/mol/K
    return {
        "n": int(len(g)), "n_temp": int(np.unique(T).size),
        "T_min": float(T.min()), "T_max": float(T.max()),
        "span_K": float(T.max() - T.min()),
        "mean_T": float(T.mean()),
        "slope": float(lr.slope), "intercept": float(lr.intercept),
        "r2": float(lr.rvalue ** 2), "p_value": float(lr.pvalue),
        "stderr_slope": float(lr.stderr),
        "dH_kJ_mol": float(dH / 1000.0),
        "dS_J_molK": float(dS),
    }


def main():
    df = load()
    rows = []
    # group by pair + solvent family; also strict same-source grouping
    for grp_keys, tag in [(("canonical_pair_id", "solvent_family"), "loose"),
                          (("canonical_pair_id", "solvent_family",
                            "source_dataset"), "strict")]:
        for key, g in df.groupby(list(grp_keys)):
            if g.temperature_K.nunique() < MIN_T:
                continue
            if g.temperature_K.max() - g.temperature_K.min() < MIN_SPAN:
                continue
            for target in ["ln_r1", "ln_r2"]:
                s = fit_series(g, target)
                if s is None:
                    continue
                rows.append({"grouping": tag,
                             "pair_id": key[0],
                             "solvent_family": key[1],
                             "source_dataset": key[2] if tag == "strict" else "",
                             "target": target,
                             "m1": g.monomer_1_name.iloc[0],
                             "m2": g.monomer_2_name.iloc[0],
                             **s})
    series = pd.DataFrame(rows)
    series.to_csv(P2 / "temperature_series_fits.csv", index=False)

    # --- summary ----------------------------------------------------------
    loose = series[series.grouping == "loose"]
    strict = series[series.grouping == "strict"]
    sig = loose[loose.p_value < 0.05]
    summary = {
        "n_loose_series": int(len(loose)),
        "n_loose_pairs": int(loose.pair_id.nunique()),
        "n_strict_series": int(len(strict)),
        "n_strict_pairs": int(strict.pair_id.nunique()),
        "n_significant_p05": int(len(sig)),
        "pct_significant_p05": round(100 * len(sig) / max(len(loose), 1), 1),
        "median_r2": round(float(loose.r2.median()), 3) if len(loose) else None,
        "median_span_K": round(float(loose.span_K.median()), 1) if len(loose) else None,
        "median_dH_kJ_mol": round(float(loose.dH_kJ_mol.median()), 2) if len(loose) else None,
        "IQR_dH_kJ_mol": [round(float(loose.dH_kJ_mol.quantile(q)), 2)
                          for q in (0.25, 0.75)] if len(loose) else None,
        "median_dS_J_molK": round(float(loose.dS_J_molK.median()), 2) if len(loose) else None,
    }

    # enthalpy-entropy compensation: dH vs dS (expect T_comp ~ 300 K)
    if len(loose) >= 5:
        lr = stats.linregress(loose.dS_J_molK, loose.dH_kJ_mol)
        summary["compensation_slope_K"] = round(float(lr.slope * 1000), 1)
        summary["compensation_r2"] = round(float(lr.rvalue ** 2), 3)

    (P2 / "p2_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("\nGo/No-Go: strict multi-T series =", summary["n_strict_series"],
          "(>=150 strong, 50-150 conditional, <50 no-go)")


if __name__ == "__main__":
    main()
