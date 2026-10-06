"""Test feed-vs-polymer composition assumptions with pair-grouped OOF CV.

The experimental Tg tables report incorporated polymer composition. The
original sequence descriptors used that composition as if it were monomer
feed. Here we also invert Mayo-Lewis to estimate the instantaneous feed that
would produce each observed incorporated composition, then compare both
feature constructions under the same leave-one-monomer-pair-out folds.

This inversion is a sensitivity analysis only: real samples may integrate
composition drift over conversion, so a single inferred feed is not a
conversion-resolved sequence reconstruction.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from sklearn.base import clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from p4_sequence import mayo_lewis, sequence_descriptors

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projects" / "P4_sequence_tg"
INPUT = PROJECT / "real_copolymer_tg_dataset.parquet"
OUTPUT_JSON = PROJECT / "p4_supplementary_feed_sensitivity.json"
OUTPUT_CSV = PROJECT / "p4_supplementary_feed_sensitivity.csv"

SEQ_FEATS = ["blockiness", "alternation", "seq_entropy", "dyad_11",
             "dyad_12", "dyad_21", "dyad_22", "P11", "P12", "P21", "P22"]
RR_FEATS = ["rA", "rB", "rA_times_rB", "ln_rA", "ln_rB"]
COMP_FEATS = ["fA", "TgA_K", "TgB_K", "deltaTg_K"]


def invert_mayo_lewis(r1: float, r2: float, polymer_fraction: float) -> float:
    target = float(polymer_fraction)
    if target <= 0.0:
        return 0.0
    if target >= 1.0:
        return 1.0

    def residual(feed: float) -> float:
        return float(np.asarray(mayo_lewis(r1, r2, feed)["F1"])) - target

    return float(brentq(residual, 1e-10, 1 - 1e-10, xtol=1e-12))


def add_inverse_feed_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    feeds = [invert_mayo_lewis(a, b, f)
             for a, b, f in zip(out.rA, out.rB, out.fA)]
    out["fA_feed_inferred"] = feeds
    desc = sequence_descriptors(out.rA.to_numpy(), out.rB.to_numpy(),
                                out.fA_feed_inferred.to_numpy())
    for name in SEQ_FEATS:
        out[f"feed_{name}"] = desc[name]
    return out


def make_features(df: pd.DataFrame, feature_set: str) -> np.ndarray:
    cols = list(COMP_FEATS)
    if feature_set == "comp_rr":
        cols += RR_FEATS
    elif feature_set == "comp_seq_polymer_composition":
        cols += SEQ_FEATS
    elif feature_set == "comp_seq_inferred_feed":
        cols += [f"feed_{c}" for c in SEQ_FEATS]
    elif feature_set != "comp":
        raise ValueError(feature_set)
    return np.nan_to_num(df[cols].to_numpy(np.float32), nan=0.0,
                         posinf=0.0, neginf=0.0)


def oof_predictions(df: pd.DataFrame, feature_set: str, model_name: str) -> np.ndarray:
    X = make_features(df, feature_set)
    y_resid = (df.Tg_exp_K - df.Tg_Fox_w_K).to_numpy(float)
    groups = df.canonical_pair_id.to_numpy()
    if model_name == "ridge":
        model = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
    elif model_name == "rf":
        model = RandomForestRegressor(n_estimators=400, random_state=0,
                                      min_samples_leaf=2, n_jobs=-1)
    else:
        raise ValueError(model_name)
    oof = np.full(len(df), np.nan)
    for train, test in LeaveOneGroupOut().split(X, y_resid, groups):
        fitted = clone(model)
        fitted.fit(X[train], y_resid[train])
        oof[test] = fitted.predict(X[test])
    return oof + df.Tg_Fox_w_K.to_numpy(float)


def score(y: np.ndarray, pred: np.ndarray) -> dict:
    err = y - pred
    denom = np.sum((y - y.mean()) ** 2)
    return {"MAE_K": float(np.mean(np.abs(err))),
            "RMSE_K": float(np.sqrt(np.mean(err ** 2))),
            "R2": float(1 - np.sum(err ** 2) / denom)}


def pair_bootstrap(df: pd.DataFrame, pred_map: dict[str, np.ndarray],
                   n_boot: int = 5000, seed: int = 71) -> dict:
    rng = np.random.default_rng(seed)
    pair_ids = df.canonical_pair_id.drop_duplicates().to_numpy()
    rows = {pid: np.flatnonzero(df.canonical_pair_id.to_numpy() == pid)
            for pid in pair_ids}
    y = df.Tg_exp_K.to_numpy(float)
    draws = {name: [] for name in pred_map}
    deltas = {name: [] for name in pred_map if name != "Fox_weight"}
    for _ in range(n_boot):
        sampled = rng.choice(pair_ids, size=len(pair_ids), replace=True)
        idx = np.concatenate([rows[pid] for pid in sampled])
        fox_mae = np.mean(np.abs(y[idx] - pred_map["Fox_weight"][idx]))
        for name, pred in pred_map.items():
            mae = np.mean(np.abs(y[idx] - pred[idx]))
            draws[name].append(float(mae))
            if name != "Fox_weight":
                deltas[name].append(float(mae - fox_mae))
    result = {}
    for name, values in draws.items():
        result[name] = {"MAE_K_CI95": np.quantile(values, [0.025, 0.975]).tolist()}
        if name in deltas:
            result[name]["MAE_minus_Fox_K_CI95"] = np.quantile(
                deltas[name], [0.025, 0.975]).tolist()
    return result


def main() -> None:
    df = pd.read_parquet(INPUT).reset_index(drop=True)
    df = add_inverse_feed_features(df)
    specs = {
        "Fox_weight": None,
        "Ridge_comp": ("comp", "ridge"),
        "Ridge_comp_rr": ("comp_rr", "ridge"),
        "Ridge_comp_seq_polymer_composition": ("comp_seq_polymer_composition", "ridge"),
        "Ridge_comp_seq_inferred_feed": ("comp_seq_inferred_feed", "ridge"),
        "RF_comp": ("comp", "rf"),
        "RF_comp_rr": ("comp_rr", "rf"),
        "RF_comp_seq_polymer_composition": ("comp_seq_polymer_composition", "rf"),
        "RF_comp_seq_inferred_feed": ("comp_seq_inferred_feed", "rf"),
    }
    y = df.Tg_exp_K.to_numpy(float)
    predictions = {"Fox_weight": df.Tg_Fox_w_K.to_numpy(float)}
    records = []
    for name, spec in specs.items():
        pred = predictions["Fox_weight"] if spec is None else oof_predictions(df, *spec)
        predictions[name] = pred
        records.append({"model": name, **score(y, pred)})

    out = {
        "n_rows": int(len(df)),
        "n_pairs": int(df.canonical_pair_id.nunique()),
        "cv": "leave-one-canonical-monomer-pair-out",
        "bootstrap": "5000 paired cluster resamples of monomer-pair groups",
        "composition_note": "observed incorporated composition remains a composition feature; inferred feed is used only for sequence descriptors",
        "inferred_feed_summary": {
            "min": float(df.fA_feed_inferred.min()),
            "median": float(df.fA_feed_inferred.median()),
            "max": float(df.fA_feed_inferred.max()),
            "median_abs_shift_from_polymer_mole_fraction": float(
                np.median(np.abs(df.fA_feed_inferred - df.fA)))
        },
        "metrics": records,
        "pair_cluster_bootstrap": pair_bootstrap(df, predictions),
        "limitations": [
            "single inferred feed is an approximation for samples that may have composition drift",
            "leave-one-pair-out can still share monomers between train and test pairs",
            "the 40-pair dataset is dominated by near-random systems and has limited power for strongly nonrandom pairs",
            "fixed model settings are inherited from the existing analysis rather than nested-hyperparameter tuned",
        ]
    }
    (PROJECT / "real_copolymer_tg_dataset_with_inferred_feed.parquet").parent.mkdir(
        parents=True, exist_ok=True)
    df.to_parquet(PROJECT / "real_copolymer_tg_dataset_with_inferred_feed.parquet", index=False)
    OUTPUT_JSON.write_text(json.dumps(out, indent=2))
    pd.DataFrame(records).to_csv(OUTPUT_CSV, index=False)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
