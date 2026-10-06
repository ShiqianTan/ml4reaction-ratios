"""P4 (real data): RR -> Mayo-Lewis sequence -> copolymer Tg delta-learning.

Uses REAL copolymer Tg datasets (radical copolymers with known composition):
  * tgvalue-0416.xlsx  'A-rand-B-Co-Tg'  (acrylates/methacrylates/St/AN/VC)
  * cleaned-2006-styrenic-co48.xlsx      (styrene + comonomer)
paired with REAL reactivity ratios from RRMaster v1.

Chain:  (r1, r2, feed mol frac f1) -> Mayo-Lewis sequence descriptors
        -> delta-learning of the Fox-mixing-rule residual.

Run: python3 src/p4_real_sequence_tg.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import STANDARDIZED, PROJECTS, R_GAS, inchikey
from p4_sequence import sequence_descriptors

P4 = PROJECTS / "P4_sequence_tg"
DL_DIR = "/Users/shiqian/Documents/project/polymer/DeltaLearningCopolymerTg"


# --------------------------------------------------------------------------- #
# reactivity-ratio lookup
# --------------------------------------------------------------------------- #

def build_rr_lookup():
    rm = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    rm = rm[rm.is_primary_measurement & (rm.r1 > 0) & (rm.r2 > 0)
            & rm.canonical_pair_id.notna()]
    g = rm.groupby("canonical_pair_id").agg(
        r1=("r1", "median"), r2=("r2", "median"), n=("r1", "size"))
    return g


def oriented_rr(ik_a, ik_b, lookup):
    """Return (r_A, r_B) so r_A=k_AA/k_AB, r_B=k_BB/k_BA."""
    if ik_a is None or ik_b is None:
        return np.nan, np.nan
    if ik_a == ik_b:
        return 1.0, 1.0
    pid = "__".join(sorted([ik_a, ik_b]))
    if pid not in lookup.index:
        return np.nan, np.nan
    r1, r2 = lookup.loc[pid, "r1"], lookup.loc[pid, "r2"]
    if ik_a <= ik_b:  # canonical order matches (m1, m2)
        return r1, r2
    return r2, r1


# --------------------------------------------------------------------------- #
# load real copolymer Tg datasets (radical copolymers)
# --------------------------------------------------------------------------- #

def load_tgvalue():
    df = pd.read_excel(f"{DL_DIR}/data/tgvalue-0416.xlsx",
                       sheet_name="A-rand-B-Co-Tg")
    df = df[df["A_random_B_Tg (Exp)"].notna()].copy()
    df["Tg_exp_K"] = df["A_random_B_Tg (Exp)"] + 273.15
    df["TgA_K"] = df["Polymer_A_Tg(°C)"] + 273.15
    df["TgB_K"] = df["Polymer_B_Tg(°C)"] + 273.15
    df["wA"] = df["Weight_Percent_A"]
    df["source"] = "tgvalue"
    return df.rename(columns={"Monomer_A_NAME": "nameA", "Monomer_B_NAME": "nameB",
                              "Monomer_A_SMILES": "smiA", "Monomer_B_SMILES": "smiB"})


def load_styrenic():
    df = pd.read_excel(f"{DL_DIR}/data/cleaned-2006-styrenic-co48.xlsx",
                       sheet_name="共聚物参数表", header=1)
    df = df[df["Tₛ (K)exp"].notna()].copy()
    df["Tg_exp_K"] = df["Tₛ (K)exp"]
    df["TgA_K"] = df["Polymer_A_Tg(°C)"] + 273.15
    df["TgB_K"] = df["Polymer_B_Tg(°C)"] + 273.15
    df["wA"] = df["Weight_Percent_A"]
    df["source"] = "styrenic"
    return df.rename(columns={"Monomer_A_NAME": "nameA", "Monomer_B_NAME": "nameB",
                              "Monomer_A_SMILES": "smiA", "Monomer_B_SMILES": "smiB"})


def mol_wt(smiles):
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    return Descriptors.MolWt(m) if m else np.nan


def build_dataset():
    lookup = build_rr_lookup()
    tg = load_tgvalue()
    st = load_styrenic()
    cols = ["source", "nameA", "nameB", "smiA", "smiB", "wA", "Tg_exp_K",
            "TgA_K", "TgB_K"]
    df = pd.concat([tg[cols], st[cols]], ignore_index=True)

    df["ikA"] = [inchikey(s) for s in df.smiA]
    df["ikB"] = [inchikey(s) for s in df.smiB]
    df["mwA"] = [mol_wt(s) for s in df.smiA]
    df["mwB"] = [mol_wt(s) for s in df.smiB]
    df["wB"] = 1 - df.wA
    # mole fraction of A in feed (from mass fraction + monomer MW)
    df["fA"] = (df.wA / df.mwA) / (df.wA / df.mwA + df.wB / df.mwB)

    rr = [oriented_rr(a, b, lookup) for a, b in zip(df.ikA, df.ikB)]
    df["rA"] = [x[0] for x in rr]
    df["rB"] = [x[1] for x in rr]
    df["canonical_pair_id"] = df.apply(
        lambda r: "__".join(sorted([r.ikA, r.ikB])) if (r.ikA and r.ikB)
        else None, axis=1)

    # keep only rows with known reactivity ratios
    matched = df[df.rA.notna() & df.rB.notna()].copy()

    # Mayo-Lewis sequence descriptors at feed composition fA
    sd = sequence_descriptors(matched.rA.values, matched.rB.values,
                              matched.fA.values)
    for k, v in sd.items():
        matched[k] = v

    # baselines (Kelvin)
    TgA, TgB = matched.TgA_K.values, matched.TgB_K.values
    wA = matched.wA.values
    matched["Tg_Fox_w_K"] = 1 / (wA / TgA + (1 - wA) / TgB)
    xA = matched.fA.values
    matched["Tg_Fox_x_K"] = 1 / (xA / TgA + (1 - xA) / TgB)
    kGT = TgA / TgB  # Gordon-Taylor (Simha-Boyer) k
    matched["Tg_GT_K"] = (wA * TgA + kGT * (1 - wA) * TgB) / (wA + kGT * (1 - wA))

    matched["resid_Fox_w_K"] = matched.Tg_exp_K - matched.Tg_Fox_w_K
    matched["deltaTg_K"] = (matched.TgA_K - matched.TgB_K).abs()
    matched["rA_times_rB"] = matched.rA * matched.rB
    matched["ln_rA"] = np.log(matched.rA)
    matched["ln_rB"] = np.log(matched.rB)
    # pair-level blockiness at equimolar feed (composition-independent)
    mid = sequence_descriptors(matched.rA.values, matched.rB.values,
                               np.full(len(matched), 0.5))
    matched["blockiness_mid"] = mid["blockiness"]
    return matched, lookup


# --------------------------------------------------------------------------- #
# delta-learning CV
# --------------------------------------------------------------------------- #

SEQ_FEATS = ["blockiness", "alternation", "seq_entropy",
             "dyad_11", "dyad_12", "dyad_21", "dyad_22",
             "P11", "P12", "P21", "P22"]
RR_FEATS = ["rA", "rB", "rA_times_rB", "ln_rA", "ln_rB"]
COMP_FEATS = ["fA", "TgA_K", "TgB_K", "deltaTg_K"]


def _ecfp(smiles, nbits=2048):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    arr = np.zeros(nbits, dtype=np.int8)
    if m:
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits)
        Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def _ecfp_block(df):
    eA = np.vstack([_ecfp(s) for s in df.smiA]).astype(np.float32)
    eB = np.vstack([_ecfp(s) for s in df.smiB]).astype(np.float32)
    return np.hstack([eA, eB])


def run_cv(df, baseline_col, feat_keys):
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestRegressor
    df = df.reset_index(drop=True)
    y_resid = (df.Tg_exp_K - df[baseline_col]).values
    groups = df.canonical_pair_id.values

    blocks = []
    if "comp" in feat_keys:
        blocks.append(df[COMP_FEATS].to_numpy(np.float32))
    if "rr" in feat_keys:
        blocks.append(df[RR_FEATS].to_numpy(np.float32))
    if "seq" in feat_keys:
        blocks.append(df[SEQ_FEATS].to_numpy(np.float32))
    if "ecfp" in feat_keys:
        blocks.append(_ecfp_block(df))
    X = np.nan_to_num(np.hstack(blocks).astype(np.float32), nan=0.0)

    logo = LeaveOneGroupOut()
    results = {}
    for mname, model in [("ridge", make_pipeline(StandardScaler(),
                                                 Ridge(alpha=1.0))),
                         ("rf", RandomForestRegressor(n_estimators=400,
                                                      random_state=0,
                                                      min_samples_leaf=2,
                                                      n_jobs=-1))]:
        from sklearn.base import clone
        oof = np.full(len(df), np.nan)
        for tr, te in logo.split(X, groups=groups):
            m = clone(model)
            m.fit(X[tr], y_resid[tr])
            oof[te] = m.predict(X[te])
        pred_Tg = oof + df[baseline_col].values
        y = df.Tg_exp_K.values
        mask = ~np.isnan(pred_Tg)
        mae = np.mean(np.abs(y[mask] - pred_Tg[mask]))
        rmse = np.sqrt(np.mean((y[mask] - pred_Tg[mask]) ** 2))
        ss = 1 - np.sum((y[mask] - pred_Tg[mask]) ** 2) / np.sum((y[mask] - y[mask].mean()) ** 2)
        results[mname] = {"MAE_K": float(mae), "RMSE_K": float(rmse),
                          "R2": float(ss)}
    return results


def main():
    df, lookup = build_dataset()
    n_pairs = df.canonical_pair_id.nunique()
    print(f"matched real copolymer Tg rows: {len(df)} across {n_pairs} pairs "
          f"(tgvalue + styrenic)")

    # baseline metrics
    def base_metrics(col):
        y = df.Tg_exp_K.values
        p = df[col].values
        mae = np.mean(np.abs(y - p))
        rmse = np.sqrt(np.mean((y - p) ** 2))
        r2 = 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
        return {"MAE_K": float(mae), "RMSE_K": float(rmse), "R2": float(r2)}

    baselines = {"Fox_weight": base_metrics("Tg_Fox_w_K"),
                 "Fox_mole": base_metrics("Tg_Fox_x_K"),
                 "GordonTaylor": base_metrics("Tg_GT_K")}

    # sequence signal on the Fox residual
    seq_signal = {}
    for f in SEQ_FEATS + ["rA_times_rB", "blockiness_mid"]:
        corr = np.corrcoef(df[f], df.resid_Fox_w_K)[0, 1]
        seq_signal[f] = round(float(corr), 3)

    # pair-level physical test: does alternation tendency (rA*rB) explain
    # the mean per-pair Fox residual?
    pair = df.groupby("canonical_pair_id").agg(
        mean_resid=("resid_Fox_w_K", "mean"),
        rA_times_rB=("rA_times_rB", "first"),
        blockiness_mid=("blockiness_mid", "first"),
        n=("resid_Fox_w_K", "size")).reset_index()
    pair_signal = {
        "corr_rArB_vs_mean_resid": round(
            float(np.corrcoef(pair.rA_times_rB, pair.mean_resid)[0, 1]), 3),
        "corr_blockiness_mid_vs_mean_resid": round(
            float(np.corrcoef(pair.blockiness_mid, pair.mean_resid)[0, 1]), 3),
        "n_pairs": int(len(pair)),
    }

    configs = {
        "delta_comp": ["comp"],
        "delta_comp+rr": ["comp", "rr"],
        "delta_comp+seq": ["comp", "seq"],
        "delta_comp+rr+ecfp": ["comp", "rr", "ecfp"],
        "delta_comp+seq+ecfp": ["comp", "seq", "ecfp"],
    }
    cv = {k: run_cv(df, "Tg_Fox_w_K", v) for k, v in configs.items()}

    out = {"n_rows": int(len(df)), "n_pairs": int(n_pairs),
           "n_tgvalue": int((df.source == "tgvalue").sum()),
           "n_styrenic": int((df.source == "styrenic").sum()),
           "baselines": baselines,
           "seq_corr_with_Fox_residual": seq_signal,
           "pair_level_alternation_signal": pair_signal,
           "delta_learning": cv}

    print(json.dumps(out, indent=2))
    (P4 / "p4_real_sequence_tg.json").write_text(json.dumps(out, indent=2))
    df.to_parquet(P4 / "real_copolymer_tg_dataset.parquet", index=False)
    print(f"\nsaved {P4/'real_copolymer_tg_dataset.parquet'} and "
          f"p4_real_sequence_tg.json")


if __name__ == "__main__":
    main()
