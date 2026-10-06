"""Reproduce a small literature benchmark: SAN NMR dyads, RR and Tg.

Inputs were transcribed from Tables 1 and 4 of Sanghvi et al. (2000).
The leave-one-out interaction coefficient is fitted only on the remaining
observed Tg rows; all reported values are exploratory because n=5 and comes
from one chemistry/source.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from p4_sequence import sequence_descriptors

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projects" / "P4_sequence_tg"
INPUT = PROJECT / "san_microemulsion_matched_case.csv"


def main():
    d = pd.read_csv(INPUT)
    d["styrene_locus_mole_fraction"] = d.locus_M0S_over_M0A / (1 + d.locus_M0S_over_M0A)
    rr = sequence_descriptors(d.r_styrene.to_numpy(float),
                              d.r_acrylonitrile.to_numpy(float),
                              d.styrene_locus_mole_fraction.to_numpy(float))
    d["FSS_RR"] = np.asarray(rr["dyad_11"])
    d["FSA_RR"] = np.asarray(rr["dyad_12"]) + np.asarray(rr["dyad_21"])
    d["FAA_RR"] = np.asarray(rr["dyad_22"])
    d["FSA_RR_minus_NMR"] = d.FSA_RR - d.FSA_NMR
    d["RR_dyad_abs_error"] = d.FSA_RR_minus_NMR.abs()

    # Fox baseline with repeat-unit mass fractions derived from measured
    # copolymer molar composition; pure Tg endpoints are reported in article.
    mw_st, mw_an = 104.15, 53.06
    w_st = (d.styrene_copolymer_mole_fraction * mw_st /
            (d.styrene_copolymer_mole_fraction * mw_st
             + (1-d.styrene_copolymer_mole_fraction) * mw_an))
    d["Fox_Tg_C"] = 1 / (w_st / (105 + 273.15)
                         + (1-w_st) / (128 + 273.15)) - 273.15

    observed = d.Tg_experimental_C.notna().to_numpy()
    o = d.loc[observed].copy().reset_index(drop=True)
    def loo_dyad_tg(prefix: str):
        fss = o[f"FSS_{prefix}"].to_numpy(float)
        fsa = o[f"FSA_{prefix}"].to_numpy(float)
        faa = o[f"FAA_{prefix}"].to_numpy(float)
        y = o.Tg_experimental_C.to_numpy(float)
        base = fss * 105 + faa * 128
        residual = y - base
        pred = np.full(len(o), np.nan)
        for i in range(len(o)):
            train = np.arange(len(o)) != i
            beta = np.sum(fsa[train] * residual[train]) / np.sum(fsa[train] ** 2)
            pred[i] = base[i] + fsa[i] * beta
        return pred

    o["NMR_dyad_LOO_Tg_C"] = loo_dyad_tg("NMR")
    o["RR_dyad_LOO_Tg_C"] = loo_dyad_tg("RR")
    fox = o.Fox_Tg_C.to_numpy(float)
    y = o.Tg_experimental_C.to_numpy(float)
    metrics = {}
    for name, pred in {
        "Fox": fox,
        "published_Barton_resubstitution": o.Tg_Barton_published_C.to_numpy(float),
        "measured_NMR_dyad_LOO": o.NMR_dyad_LOO_Tg_C.to_numpy(float),
        "RR_derived_dyad_LOO": o.RR_dyad_LOO_Tg_C.to_numpy(float),
    }.items():
        err = pred - y
        metrics[name] = {"n": int(len(y)), "MAE_C": float(np.mean(np.abs(err))),
                         "RMSE_C": float(np.sqrt(np.mean(err**2)))}
    d.to_csv(PROJECT / "san_microemulsion_matched_case_results.csv", index=False)
    o.to_csv(PROJECT / "san_microemulsion_tg_loo_predictions.csv", index=False)
    out = {
        "n_rows_transcribed": int(len(d)),
        "n_rows_with_experimental_tg": int(len(o)),
        "RR_predicted_heterodyad_MAE_all_rows": float(d.RR_dyad_abs_error.mean()),
        "RR_predicted_heterodyad_MAE_experimental_Tg_rows": float(o.RR_dyad_abs_error.mean()),
        "Tg_metrics": metrics,
        "conditions": "Styrene/acrylonitrile in anionic oil-in-water microemulsion; 70 C; 5% SAN monomer, 14% SDS, 81% water; KPS 0.37 mM.",
        "cautions": [
            "only five rows have experimental Tg and all are from one monomer pair and one study",
            "RR are source-reported system-level values; the held-out Tg is not held out from RR estimation",
            "the locus-feed fraction is reconstructed from the reported M0S/M0A ratio",
            "the heterodyad-based Tg interaction coefficient is refit within the five-row set; this is not cross-pair validation",
            "published Barton values were calculated from measured NMR dyads and should not be treated as RR-only predictions",
        ],
    }
    (PROJECT / "p4_san_matched_case.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
