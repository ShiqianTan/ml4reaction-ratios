"""P4: end-to-end RR -> sequence -> Tg validation on real sequence-Tg data.

Uses the curated experimental sequence-Tg dataset from the sibling repo
``/Users/shiqian/Documents/project/polymer/SequenceTg`` (Johnston 1973/1976,
Barton 1970, plus NMR-dyad radical systems), which supplies the strongly
alternating / sequence-controlled radical systems that the DeltaLearningCopolymerTg
data lacked.

Two parts:
  A. RR -> sequence: validate the Mayo-Lewis engine against NMR-measured dyads.
  B. RR -> sequence -> Tg: Barton dyad-additivity model (dyads from Mayo-Lewis,
     Bernoulli ideal-random, or measured) vs the composition-only Fox baseline.

Run: python3 src/p4_sequence_tg_validation.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import PROJECTS
from p4_sequence import sequence_descriptors

P4 = PROJECTS / "P4_sequence_tg"
SEQTG = "/Users/shiqian/Documents/project/polymer/SequenceTg"

K = 273.15


def fox_w(w_A, tA, tB):
    """Fox equation on weight fractions (temperatures in K)."""
    return 1.0 / (w_A / tA + (1 - w_A) / tB)


def barton(dyads, tAA, tBB, tAB):
    """Barton dyad-additivity: Tg = n_AA Tg_AA + n_BB Tg_BB + n_AB Tg_AB (K)."""
    return dyads[0] * tAA + dyads[1] * tBB + (dyads[2] + dyads[3]) * tAB


def bernoulli_dyads(fA):
    fB = 1 - fA
    return (fA ** 2, fB ** 2, fA * fB, fA * fB)


def ml_dyads(r1, r2, fA):
    s = sequence_descriptors(r1, r2, fA)
    return (s["dyad_11"], s["dyad_22"], s["dyad_12"], s["dyad_21"])


def main():
    out = {}

    # ------------------------------------------------------------------ #
    # Part A: RR -> sequence (Mayo-Lewis vs NMR dyads)
    # r1/r2 oriented to (monomer_A, monomer_B); NMR dyads = (AA, BB, AB)
    # ------------------------------------------------------------------ #
    systems = [
        dict(name="AN-EA", r1=0.92, r2=1.20, fA=0.50, nmr=(0.33, 0.19, 0.48)),
        dict(name="BVE-CEVE", r1=0.88, r2=0.92, fA=0.50, nmr=(0.28, 0.24, 0.48)),
        dict(name="NVP-CEVE", r1=8.12, r2=0.11, fA=0.50, nmr=(0.81, 0.01, 0.18)),
        dict(name="S-Acrylamide", r1=1.04, r2=0.98, fA=0.713, nmr=(0.50, 0.10, 0.40)),
        dict(name="MMA-S(random)", r1=0.46, r2=0.52, fA=0.50, nmr=(0.25, 0.25, 0.50)),
        dict(name="VDF-HFP", r1=2.45, r2=0.00, fA=0.89, nmr=(0.78, 0.11, 0.11)),
    ]
    rowsA = []
    for s in systems:
        d = ml_dyads(s["r1"], s["r2"], s["fA"])
        pred = (d[0], d[1], d[2] + d[3])
        err = np.array(pred) - np.array(s["nmr"])
        rowsA.append(dict(system=s["name"], r1=s["r1"], r2=s["r2"], fA=s["fA"],
                          nmr_AA=s["nmr"][0], nmr_BB=s["nmr"][1],
                          nmr_AB=s["nmr"][2],
                          pred_AA=round(pred[0], 3), pred_BB=round(pred[1], 3),
                          pred_AB=round(pred[2], 3),
                          err_AA=round(err[0], 3), err_BB=round(err[1], 3),
                          err_AB=round(err[2], 3)))
    A = pd.DataFrame(rowsA)
    A.to_csv(P4 / "rr_sequence_dyad_validation.csv", index=False)
    mae_dyad = float(np.mean(np.abs(A[["err_AA", "err_BB", "err_AB"]]).to_numpy()))
    out["partA"] = {"n_systems": len(A), "dyad_MAE": round(mae_dyad, 3),
                    "hetero_dyad_MAE": round(float(A.err_AB.abs().mean()), 3)}
    print("=== Part A: RR -> sequence (Mayo-Lewis vs NMR dyads) ===")
    print(A[["system", "nmr_AB", "pred_AB", "err_AB"]].round(3).to_string(index=False))
    print(f"dyad MAE = {mae_dyad:.3f}\n")

    # ------------------------------------------------------------------ #
    # Part B: RR -> sequence -> Tg on Johnston 1973/1976 systems
    # ------------------------------------------------------------------ #
    j = pd.read_csv(f"{SEQTG}/data/johnston_sequence_tg_master.csv")
    # reactivity ratios (monomer_A, monomer_B); None -> use Bernoulli
    rr = {"aMS-AN": (0.17, 0.088),      # ZnCl2-complexed (alternating)
          "MMA-VCl": (11.2, 0.044),
          "BMA-VCl": None,
          "S-AN": (0.40, 0.07),         # RRMaster median
          "AN-MMA": (0.14, 1.30)}       # RRMaster median

    rowsB = []
    for _, r in j.iterrows():
        fA = r.mol_A / 100.0
        wA = r.mol_A * r.M_A / (r.mol_A * r.M_A + (100 - r.mol_A) * r.M_B)
        tAA, tBB, tAB = r.Tg_homo_A_C + K, r.Tg_homo_B_C + K, r.Tg_AB_C + K
        t_exp = r.Tg_C + K
        fox = fox_w(wA, tAA, tBB)
        # dyad model: alternating sample -> pure AB; random -> Bernoulli (or ML)
        if r.seq_type == "alternating":
            pred = tAB
            mode = "alternating"
        else:
            rr_sys = rr.get(r.system)
            if rr_sys is not None:
                dy = ml_dyads(rr_sys[0], rr_sys[1], fA)
                pred = barton(dy, tAA, tBB, tAB)
                mode = "mlewis"
            else:
                pred = barton(bernoulli_dyads(fA), tAA, tBB, tAB)
                mode = "bernoulli"
        rowsB.append(dict(system=r.system, seq_type=r.seq_type, mol_A=r.mol_A,
                          Tg_exp_C=r.Tg_C, fox_C=fox - K, dyad_C=pred - K,
                          mode=mode,
                          err_fox=round((t_exp - fox), 1),
                          err_dyad=round((t_exp - pred), 1)))
    B = pd.DataFrame(rowsB)
    B.to_csv(P4 / "rr_sequence_tg_validation.csv", index=False)
    mae_fox = B.err_fox.abs().mean()
    mae_dyad = B.err_dyad.abs().mean()
    r2_fox = 1 - ((B.Tg_exp_C - B.fox_C) ** 2).sum() / ((B.Tg_exp_C - B.Tg_exp_C.mean()) ** 2).sum()
    r2_dyad = 1 - ((B.Tg_exp_C - B.dyad_C) ** 2).sum() / ((B.Tg_exp_C - B.Tg_exp_C.mean()) ** 2).sum()
    out["partB"] = {"n_samples": len(B),
                    "fox_MAE_K": round(float(mae_fox), 1),
                    "dyad_MAE_K": round(float(mae_dyad), 1),
                    "fox_R2": round(float(r2_fox), 3),
                    "dyad_R2": round(float(r2_dyad), 3)}
    print("=== Part B: RR -> sequence -> Tg (Barton dyad vs Fox) ===")
    for sys, g in B.groupby("system"):
        print(f"  {sys:9s} n={len(g):2d}  Fox MAE={g.err_fox.abs().mean():5.1f}K  "
              f"dyad MAE={g.err_dyad.abs().mean():5.1f}K")
    print(f"  OVERALL  Fox MAE={mae_fox:.1f}K (R2={r2_fox:.3f})  "
          f"dyad MAE={mae_dyad:.1f}K (R2={r2_dyad:.3f})\n")

    # ------------------------------------------------------------------ #
    # Part C: Barton sign rule across systems (Tg_AB vs homopolymer mean)
    # ------------------------------------------------------------------ #
    sub = j.drop_duplicates("system")
    rowsC = []
    for _, r in sub.iterrows():
        avg = (r.Tg_homo_A_C + r.Tg_homo_B_C) / 2.0
        rowsC.append(dict(system=r.system, Tg_AB_C=r.Tg_AB_C, avg_homo_C=round(avg, 1),
                          dTg_AB_minus_avg=round(r.Tg_AB_C - avg, 1),
                          shape="Min" if r.Tg_AB_C < avg else "Max"))
    C = pd.DataFrame(rowsC)
    C.to_csv(P4 / "barton_sign_rule.csv", index=False)
    out["partC"] = {"alternating_lowers_Tg": C.shape[0] and
                    int((C.dTg_AB_minus_avg < 0).sum()),
                    "alternating_raises_Tg": int((C.dTg_AB_minus_avg > 0).sum())}
    print("=== Part C: Barton sign rule (alternating dyad Tg vs homopolymer mean) ===")
    print(C.to_string(index=False))

    (P4 / "p4_sequence_tg_validation.json").write_text(json.dumps(out, indent=2))
    print(f"\nWrote outputs to {P4}/  summary: {json.dumps(out, indent=2)}")


if __name__ == "__main__":
    main()
