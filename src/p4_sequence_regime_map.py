"""Exploratory P4 regime map: sequence leverage vs Fox composition baseline.

For each Johnston sample, estimate the change in Barton Tg caused by moving
from a Bernoulli sequence at fixed incorporated composition to the sequence
specified by the table (alternating samples) or by terminal RR (random samples):

    dTg_seq = (H_seq - 2*xA*(1-xA)) *
              (Tg_AB - (Tg_AA + Tg_BB)/2)

This is a diagnostic decomposition, not held-out prediction: Tg_AB is a
system-level parameter from the same source, and some RR are source-level or
database estimates. The output is intended to identify high-leverage regimes
and expose where sequence leverage is too small to plausibly explain a large
Fox residual.

Run: python3 src/p4_sequence_regime_map.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from p4_sequence import mayo_lewis, sequence_descriptors

ROOT = Path(__file__).resolve().parents[1]
P4 = ROOT / "projects" / "P4_sequence_tg"
INPUT = Path("/Users/shiqian/Documents/project/polymer/SequenceTg/data/johnston_sequence_tg_master.csv")
OUT_CSV = P4 / "p4_sequence_regime_map.csv"
OUT_SUMMARY = P4 / "p4_sequence_regime_map.json"
K = 273.15

# Keep aligned with the existing Johnston RR->Tg audit. Values are
# source-level/model-median estimates, not fold-specific fits.
RR = {
    "aMS-AN": (0.17, 0.088),
    "MMA-VCl": (11.2, 0.044),
    "BMA-VCl": None,
    "S-AN": (0.40, 0.07),
    "AN-MMA": (0.14, 1.30),
}


def invert_mayo_lewis(r1: float, r2: float, target_f: float) -> float:
    """Bisection inversion of instantaneous Mayo-Lewis copolymer fraction."""
    lo, hi = 1e-10, 1 - 1e-10
    for _ in range(100):
        mid = (lo + hi) / 2
        f_out = float(np.asarray(mayo_lewis(r1, r2, mid)["F1"]))
        if f_out < target_f:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def fox_tg_c(row: pd.Series) -> float:
    w_a = row.mol_A * row.M_A / (
        row.mol_A * row.M_A + (100 - row.mol_A) * row.M_B
    )
    tg_a, tg_b = row.Tg_homo_A_C + K, row.Tg_homo_B_C + K
    return 1.0 / (w_a / tg_a + (1 - w_a) / tg_b) - K


def main() -> None:
    df = pd.read_csv(INPUT)
    rows: list[dict] = []
    for _, row in df.iterrows():
        x_a = float(row.mol_A / 100)
        h_random = 2 * x_a * (1 - x_a)
        rr = RR.get(row.system)

        if row.seq_type == "alternating":
            h_sequence = 1.0
            sequence_basis = "reported alternating class"
            rr_status = "not used for this row"
        elif rr is None:
            h_sequence = h_random
            sequence_basis = "Bernoulli fallback"
            rr_status = "RR unavailable; no sequence leverage inferred"
        else:
            feed_a = invert_mayo_lewis(rr[0], rr[1], x_a)
            desc = sequence_descriptors(rr[0], rr[1], feed_a)
            h_sequence = float(np.asarray(desc["dyad_12"] + desc["dyad_21"]))
            sequence_basis = "terminal RR, inverted to reported polymer composition"
            rr_status = f"system-level RR {rr[0]:g}, {rr[1]:g}"

        delta_h = h_sequence - h_random
        delta_t_ab = float(
            row.Tg_AB_C - (row.Tg_homo_A_C + row.Tg_homo_B_C) / 2
        )
        sequence_shift = delta_h * delta_t_ab
        tg_fox = fox_tg_c(row)
        rows.append({
            "system": row.system,
            "source": row.source,
            "seq_type": row.seq_type,
            "monomer_A_mol_pct": row.mol_A,
            "Tg_exp_C": row.Tg_C,
            "Tg_Fox_C": tg_fox,
            "observed_Fox_residual_K": row.Tg_C - tg_fox,
            "heterodyad_random": h_random,
            "heterodyad_sequence_model": h_sequence,
            "heterodyad_shift_vs_random": delta_h,
            "Tg_AB_minus_homo_mean_K": delta_t_ab,
            "Barton_sequence_shift_vs_random_K": sequence_shift,
            "sequence_basis": sequence_basis,
            "RR_status": rr_status,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    summary = []
    for system, group in out.groupby("system", sort=True):
        nonzero = group[group.sequence_basis != "Bernoulli fallback"]
        summary.append({
            "system": system,
            "n_samples": int(len(group)),
            "n_rows_with_nontrivial_sequence_basis": int(len(nonzero)),
            "Tg_AB_minus_homo_mean_K": float(group.Tg_AB_minus_homo_mean_K.iloc[0]),
            "median_abs_heterodyad_shift": float(group.heterodyad_shift_vs_random.abs().median()),
            "median_abs_Barton_sequence_shift_K": float(group.Barton_sequence_shift_vs_random_K.abs().median()),
            "range_Barton_sequence_shift_K": [
                float(group.Barton_sequence_shift_vs_random_K.min()),
                float(group.Barton_sequence_shift_vs_random_K.max()),
            ],
            "Fox_MAE_K": float(group.observed_Fox_residual_K.abs().mean()),
            "interpretation": (
                "diagnostic only; Tg_AB is not held out and RR are source-level estimates"
            ),
        })

    result = {
        "analysis": "fixed-composition sequence leverage under Barton dyad additivity",
        "n_samples": int(len(out)),
        "n_systems": int(out.system.nunique()),
        "equation": "(H_sequence - 2*xA*(1-xA)) * (Tg_AB - (Tg_AA+Tg_BB)/2)",
        "composition_note": "reported incorporated polymer composition; RR feed inverted to match it for this diagnostic",
        "limitations": [
            "Tg_AB is system-level and not held out from the samples being summarized",
            "RR values are system-level values aligned to the existing validation script, not uncertainty distributions",
            "BMA-VCl has no RR input and therefore uses a Bernoulli fallback; its zero sequence shift is not evidence of no sequence effect",
            "Fox residual includes composition-rule and chemistry effects, so it is not a pure sequence residual",
            "computed shifts are sensitivities under Barton additivity, not independent predictions or causal estimates",
        ],
        "systems": summary,
    }
    OUT_SUMMARY.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(f"\nWrote {OUT_CSV} and {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
