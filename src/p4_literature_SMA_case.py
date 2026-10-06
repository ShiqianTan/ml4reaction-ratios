"""Reproduce the Brar & Sunita styrene/methyl acrylate literature case.

Table values are transcribed from Tables 3-4 of the 1992 article. It reports
terminal Kelen-Tudos reactivity ratios, NMR dyads, Barton Tg estimates, and
DSC Tg for three samples. Results are an exploratory single-source case study.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from p4_sequence import mayo_lewis, sequence_descriptors

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projects" / "P4_sequence_tg"
INPUT = PROJECT / "styrene_methyl_acrylate_matched_case.csv"

MW_ST, MW_MA = 104.15, 86.09
TG_ST_K, TG_MA_K = 373.15, 279.15  # reported homopolymer values: 100 and 6 C


def invert_copolymer_composition(x_st: float, r_st: float, r_ma: float) -> float:
    """Find locus/feed f_st whose terminal Mayo-Lewis F_st matches measured x."""
    lo, hi = 1e-10, 1 - 1e-10
    for _ in range(100):
        mid = (lo + hi) / 2
        F = float(mayo_lewis(r_st, r_ma, mid)["F1"])
        if F < x_st:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def add_dyads(d: pd.DataFrame, suffix: str, desc: dict) -> None:
    d[f"FSS_{suffix}"] = np.asarray(desc["dyad_11"], float)
    d[f"FSA_{suffix}"] = (np.asarray(desc["dyad_12"], float)
                           + np.asarray(desc["dyad_21"], float))
    d[f"FAA_{suffix}"] = np.asarray(desc["dyad_22"], float)


def loo_interaction(o: pd.DataFrame, suffix: str) -> np.ndarray:
    """Leave one Tg out; fit only one Tg_SM value with fixed homo endpoints."""
    fss = o[f"FSS_{suffix}"].to_numpy(float)
    fsa = o[f"FSA_{suffix}"].to_numpy(float)
    faa = o[f"FAA_{suffix}"].to_numpy(float)
    y = o.Tg_experimental_K.to_numpy(float)
    base = fss * TG_ST_K + faa * TG_MA_K
    residual = y - base
    pred = np.full(len(o), np.nan)
    for i in range(len(o)):
        train = np.arange(len(o)) != i
        den = np.sum(fsa[train] ** 2)
        beta = np.sum(fsa[train] * residual[train]) / den
        pred[i] = base[i] + fsa[i] * beta
    return pred


def metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    err = np.asarray(pred) - np.asarray(y)
    return {"n": int(len(y)), "MAE_K": float(np.mean(np.abs(err))),
            "RMSE_K": float(np.sqrt(np.mean(err**2)))}


def main() -> None:
    d = pd.read_csv(INPUT)
    d["styrene_polymer_mole_fraction_from_NMR"] = d.FSS_NMR + 0.5 * d.FSA_NMR
    d["styrene_composition_from_NMR_dyads"] = d.styrene_polymer_mole_fraction_from_NMR

    rr_feed = sequence_descriptors(d.r_styrene.to_numpy(float),
                                   d.r_methyl_acrylate.to_numpy(float),
                                   d.styrene_feed_mole_fraction.to_numpy(float))
    add_dyads(d, "RR_at_reported_feed", rr_feed)
    d["RR_at_reported_feed_predicted_styrene_fraction"] = rr_feed["F1"]
    d["RR_at_reported_feed_heterodyad_abs_error"] = (
        d.FSA_RR_at_reported_feed - d.FSA_NMR).abs()

    x_st = d.styrene_polymer_mole_fraction_from_NMR.to_numpy(float)
    f_match = np.array([invert_copolymer_composition(x, r1, r2)
                        for x, r1, r2 in zip(x_st, d.r_styrene, d.r_methyl_acrylate)])
    d["feed_fraction_inverted_from_NMR_composition"] = f_match
    rr_fixed_x = sequence_descriptors(d.r_styrene.to_numpy(float),
                                      d.r_methyl_acrylate.to_numpy(float), f_match)
    add_dyads(d, "RR_at_measured_composition", rr_fixed_x)
    d["RR_at_measured_composition_heterodyad_abs_error"] = (
        d.FSA_RR_at_measured_composition - d.FSA_NMR).abs()

    w_st = x_st * MW_ST / (x_st * MW_ST + (1 - x_st) * MW_MA)
    d["Fox_Tg_K"] = 1 / (w_st / TG_ST_K + (1 - w_st) / TG_MA_K)

    o = d[d.Tg_experimental_K.notna()].copy().reset_index(drop=True)
    y = o.Tg_experimental_K.to_numpy(float)
    o["NMR_dyad_LOO_Tg_K"] = loo_interaction(o, "NMR")
    o["RR_reported_feed_dyad_LOO_Tg_K"] = loo_interaction(o, "RR_at_reported_feed")
    o["RR_measured_composition_dyad_LOO_Tg_K"] = loo_interaction(
        o, "RR_at_measured_composition")

    out = {
        "source": "Brar & Sunita, Indian Journal of Chemistry 31A (1992), 903-908",
        "n_transcribed": int(len(d)), "n_with_experimental_Tg": int(len(o)),
        "RR": {"r_styrene": 1.20, "r_methyl_acrylate": 0.18,
               "method": "Kelen-Tudos", "product": 0.216},
        "conditions": "bulk polymerization at 70 C, BPO 0.5 wt%, 5-10% conversion",
        "RR_heterodyad_MAE_vs_NMR_at_reported_feed": float(
            d.RR_at_reported_feed_heterodyad_abs_error.mean()),
        "RR_heterodyad_MAE_vs_NMR_composition_constrained": float(
            d.RR_at_measured_composition_heterodyad_abs_error.mean()),
        "measured_Tg_LOO_metrics": {
            "Fox_from_NMR_dyad_composition": metrics(y, o.Fox_Tg_K),
            "measured_NMR_dyad_interaction": metrics(y, o.NMR_dyad_LOO_Tg_K),
            "RR_dyad_at_reported_feed": metrics(y, o.RR_reported_feed_dyad_LOO_Tg_K),
            "RR_dyad_composition_constrained": metrics(
                y, o.RR_measured_composition_dyad_LOO_Tg_K),
        },
        "cautions": [
            "only three samples have measured DSC Tg, all from one paper and one monomer pair",
            "RR are system-level estimates from the same low-conversion study, not refit within each Tg fold",
            "measured NMR dyad and feed compositions are used as reported; NMR-derived composition follows xS=FSS+FSA/2",
            "composition-constrained RR dyads invert Mayo-Lewis F1 to match NMR-derived polymer composition; this is a sensitivity analysis, not independent feed measurement",
            "one heterodyad Tg interaction parameter is fit on two Tg observations per fold; all LOO errors are highly unstable at n=3",
            "Table 4 Barton values are included for audit but are not independent predictions because they use measured sequence data",
        ],
    }
    d.to_csv(PROJECT / "styrene_methyl_acrylate_matched_case_results.csv", index=False)
    o.to_csv(PROJECT / "styrene_methyl_acrylate_tg_loo_predictions.csv", index=False)
    (PROJECT / "p4_sma_matched_case.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
