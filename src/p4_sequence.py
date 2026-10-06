"""P4a: Mayo-Lewis sequence engine  (r1, r2, f1) -> composition + sequence.

Instantaneous copolymerisation statistics for a binary system:
  F1 = (r1 f1^2 + f1 f2) / (r1 f1^2 + 2 f1 f2 + r2 f2^2)
  P(1|1) = r1 f1 / (r1 f1 + f2)
  P(2|2) = r2 f2 / (f1 + r2 f2)
  L1 = 1/(1-P(1|1)) = 1/P(1|2), L2 = 1/P(2|1)

Sequence descriptors: dyad/triad probabilities, run lengths, blockiness,
alternation, sequence entropy.  Validated by Monte-Carlo simulation.

Run: python3 src/p4_sequence.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import STANDARDIZED, PROJECTS

P4 = PROJECTS / "P4_sequence_tg"


def mayo_lewis(r1, r2, f1):
    r1, r2, f1 = np.asarray(r1, float), np.asarray(r2, float), np.asarray(f1, float)
    f2 = 1 - f1
    denom = r1 * f1 ** 2 + 2 * f1 * f2 + r2 * f2 ** 2
    F1 = (r1 * f1 ** 2 + f1 * f2) / denom
    P11 = r1 * f1 / (r1 * f1 + f2)
    P22 = r2 * f2 / (f1 + r2 * f2)
    P12 = 1 - P11
    P21 = 1 - P22
    return {"F1": F1, "P11": P11, "P12": P12, "P21": P21, "P22": P22}


def sequence_descriptors(r1, r2, f1):
    s = mayo_lewis(r1, r2, f1)
    P11, P12, P21, P22 = s["P11"], s["P12"], s["P21"], s["P22"]
    F1 = s["F1"]
    f1 = np.asarray(f1, float)
    f2 = 1 - f1
    # dyad fractions (terminal model)
    dyad_11 = F1 * P11
    dyad_12 = F1 * P12
    dyad_21 = (1 - F1) * P21
    dyad_22 = (1 - F1) * P22
    tot = dyad_11 + dyad_12 + dyad_21 + dyad_22
    blockiness = P11 + P22          # 1 = ideal random, >1 blocky, <1 alternating
    alternation = P12 + P21         # 2 - blockiness, 1 = ideal random
    L1 = 1.0 / np.clip(P12, 1e-9, None)   # number-avg run length of 1
    L2 = 1.0 / np.clip(P21, 1e-9, None)
    # sequence (dyad) entropy, normalised
    props = np.vstack([dyad_11, dyad_12, dyad_21, dyad_22]) / tot
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -np.nansum(np.where(props > 0, props * np.log2(props), 0), axis=0)
    return {
        "F1": F1, "weight_1": F1, "P11": P11, "P12": P12, "P21": P21,
        "P22": P22, "blockiness": blockiness, "alternation": alternation,
        "L1": L1, "L2": L2, "seq_entropy": H,
        "dyad_11": dyad_11, "dyad_12": dyad_12,
        "dyad_21": dyad_21, "dyad_22": dyad_22,
    }


def _mc_validate(r1, r2, f1, n=200000, seed=0):
    """Monte-Carlo terminal-model chain to check analytic descriptors."""
    rng = np.random.default_rng(seed)
    # Markov chain on radical end state (0 == ends in monomer 1)
    state = 0
    chain = []
    for _ in range(n + 2000):
        if state == 0:
            p1 = r1 * f1 / (r1 * f1 + (1 - f1))
            state = 0 if rng.random() < p1 else 1
        else:
            p2 = r2 * (1 - f1) / (f1 + r2 * (1 - f1))
            state = 1 if rng.random() < p2 else 0
        chain.append(state)
    c = np.array(chain[2000:])
    F1_mc = float(np.mean(c == 0))   # fraction of monomer-1 units
    trans_11 = np.mean((c[:-1] == 0) & (c[1:] == 0)) / max(np.mean(c[:-1] == 0), 1e-9)
    trans_22 = np.mean((c[:-1] == 1) & (c[1:] == 1)) / max(np.mean(c[:-1] == 1), 1e-9)
    return F1_mc, trans_11, trans_22


def main():
    de = sequence_descriptors(1.0, 2.0, 0.5)
    F1_mc, t11, t22 = _mc_validate(1.0, 2.0, 0.5)
    check = {"F1_analytic": float(de["F1"]), "F1_mc": float(F1_mc),
             "P11_analytic": float(de["P11"]), "P11_mc": float(t11),
             "P22_analytic": float(de["P22"]), "P22_mc": float(t22)}
    print("Monte-Carlo validation (r1=1,r2=2,f1=0.5):")
    print(json.dumps(check, indent=2))

    # apply to real RR pairs from RRMaster
    df = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    df = df[df.ln_r1.notna() & df.ln_r2.notna()].drop_duplicates("canonical_pair_id")
    f0 = np.linspace(0.05, 0.95, 19)
    rows = []
    for _, r in df.iterrows():
        for f in f0:
            d = sequence_descriptors(r.r1, r.r2, f)
            rows.append({"canonical_pair_id": r.canonical_pair_id,
                         "m1": r.monomer_1_name, "m2": r.monomer_2_name,
                         "f1_feed": f, **d})
    seq = pd.DataFrame(rows)
    seq.to_parquet(P4 / "pair_sequence_descriptors.parquet", index=False)
    print(f"\nsequence descriptors for {df.canonical_pair_id.nunique()} pairs "
          f"x {len(f0)} feed compositions -> {len(seq)} rows")
    print(seq[["F1", "blockiness", "alternation", "L1", "L2",
               "seq_entropy"]].describe().round(3).to_string())
    (P4 / "sequence_mc_validation.json").write_text(json.dumps(check, indent=2))


if __name__ == "__main__":
    main()
