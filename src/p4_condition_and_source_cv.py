"""Condition-stratified RR -> dyad spread and source-held-out Tg CV.

This is a descriptive audit, not a causal decomposition. Categorical condition
terms are evaluated after pair fixed effects so pair chemistry is not mistaken
for within-pair condition variation. Source-held-out Tg results are reported
both with and without monomer-pair overlap across Tg sources.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestRegressor

from common import STANDARDIZED
from p4_sequence import sequence_descriptors

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projects" / "P4_sequence_tg"


def _category(s: pd.Series) -> pd.Series:
    return (s.fillna("unknown").astype(str).str.strip().str.lower()
            .replace({"": "unknown", "nan": "unknown", "na": "unknown",
                      "none": "unknown", "not reported": "unknown"}))


def _fit_method(s: pd.Series) -> pd.Series:
    x = _category(s)
    def norm(v: str) -> str:
        if v == "unknown":
            return v
        t = re.sub(r"[^a-z]+", " ", v).strip()
        for key, label in [("fineman ross", "fineman-ross"),
                           ("kelen tud", "kelen-tudos"),
                           ("extended fr", "extended-fineman-ross"),
                           ("error in variables", "errors-in-variables"),
                           ("evm", "errors-in-variables"),
                           ("nonlinear least", "nonlinear-regression"),
                           ("non linear least", "nonlinear-regression"),
                           ("mayo lewis", "mayo-lewis"),
                           ("curve fit", "nonlinear-regression")]:
            if key in t:
                return label
        return t[:60] or "unknown"
    return x.map(norm)


def _pair_fixed_effect_share(y: np.ndarray, pair: pd.Series,
                             factor: pd.Series) -> dict:
    d = pd.DataFrame({"y": y, "pair": pair.astype(str).values,
                      "factor": _category(factor).values})
    d = d.dropna(subset=["y", "pair"])
    # Missing condition labels do not constitute a physical condition level.
    d = d[d.factor != "unknown"].copy()
    # Require replicated pair-condition cells; singleton cells would turn
    # their outcome into a perfect category fit and exaggerate explained SS.
    cell_counts = d.groupby(["pair", "factor"]).size()
    keep = cell_counts[cell_counts >= 2].index
    d = d.set_index(["pair", "factor"]).loc[keep].reset_index()
    if len(d) < 3 or d.factor.nunique() < 2:
        return {"n_rows": int(len(d)), "n_strata": int(d.pair.nunique()),
                "within_pair_ss": None, "factor_explained_fraction": None}
    pair_mean = d.groupby("pair").y.transform("mean")
    total = float(np.sum((d.y.to_numpy(float) - pair_mean.to_numpy(float)) ** 2))
    cell = d.groupby(["pair", "factor"]).y.agg(["mean", "size"]).reset_index()
    pmean = d.groupby("pair").y.mean().to_dict()
    cell_ss = float(sum(row["size"] * (row["mean"] - pmean[row["pair"]]) ** 2
                        for _, row in cell.iterrows()))
    sse0, sse1 = total, max(0., total - cell_ss)
    return {"n_rows": int(len(d)), "n_strata": int(d.pair.nunique()),
            "n_factor_levels": int(d.factor.nunique()),
            "within_pair_ss": sse0,
            "factor_explained_fraction": float(max(0., (sse0 - sse1) / sse0))
            if sse0 > 0 else None,
            "singleton_pair_factor_cells_excluded": int((cell_counts < 2).sum()),
        "factor_explained_fraction_including_unknown": False}


def condition_spread() -> dict:
    tg = pd.read_parquet(PROJECT / "real_copolymer_tg_dataset.parquet")
    rr = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    rr = rr[rr.is_primary_measurement & (rr.r1 > 0) & (rr.r2 > 0)
            & rr.canonical_pair_id.isin(tg.canonical_pair_id.unique())].copy()
    # Fixed common composition removes differences in Tg sample compositions.
    ds = sequence_descriptors(rr.r1.to_numpy(float), rr.r2.to_numpy(float),
                              np.full(len(rr), .5))
    rr["heterodyad_f05"] = np.asarray(ds["dyad_12"]) + np.asarray(ds["dyad_21"])
    bins = [-np.inf, 298.5, 323.5, 348.5, np.inf]
    rr["temperature_level"] = pd.cut(rr.temperature_nominal_K, bins=bins,
                                      labels=["<=298K", "299-323K", "324-348K", ">=349K"]).astype("object")
    rr["temperature_level"] = rr.temperature_level.fillna("unknown")
    rr["source_level"] = _category(rr.source_dataset)
    ref = _category(rr.source_reference)
    doi = _category(rr.source_doi)
    rr["citation_source"] = ref.where(ref != "unknown", doi)
    rr["citation_source"] = _category(rr.citation_source)
    rr["solvent_level"] = _category(rr.solvent_family)
    rr["fitting_level"] = _fit_method(rr.fitting_method)
    fit_map = {"kelen-tudos": "Kelen-Tudos", "fineman-ross": "Fineman-Ross",
               "finemann ross": "Fineman-Ross", "errors-in-variables": "Errors-in-variables",
               "nonlinear-regression": "Nonlinear regression", "mayo-lewis": "Mayo-Lewis"}
    rr["fitting_level"] = rr.fitting_level.map(lambda x: fit_map.get(x, "unknown" if x == "unknown" else "Other"))
    rr["polymerization_level"] = _category(rr.polymerization_method)

    factors = {
        "nominal_temperature_K": "temperature_level",
        "solvent_family": "solvent_level",
        "fitting_method_grouped": "fitting_level",
        "polymerization_method": "polymerization_level",
    }
    complete = {}
    for out, col in factors.items():
        vc = rr[col].value_counts(dropna=False)
        complete[out] = {
            "known_rows": int((rr[col] != "unknown").sum()),
            "rows": int(len(rr)),
            "known_fraction": float((rr[col] != "unknown").mean()),
            "levels_including_unknown": int(rr[col].nunique()),
            "top_levels": {str(k): int(v) for k, v in vc.head(12).items()},
        }
    decomp = {}
    decomp_source_controlled = {}
    cited = rr[rr.citation_source != "unknown"].copy()
    cited["pair_source_group"] = (cited.canonical_pair_id.astype(str)
                                  + "||" + cited.citation_source.astype(str))
    for out, col in factors.items():
        decomp[out] = _pair_fixed_effect_share(rr.heterodyad_f05.to_numpy(),
                                                rr.canonical_pair_id, rr[col])
        decomp_source_controlled[out] = _pair_fixed_effect_share(
            cited.heterodyad_f05.to_numpy(), cited.pair_source_group, cited[col])
    # Cross-citation pair-level dyad spreads at common f=0.5, weighted equally
    # by citation within a pair. Corpus buckets are not treated as sources.
    med = rr[rr.citation_source != "unknown"].groupby(
        ["canonical_pair_id", "citation_source"], as_index=False).agg(
        dyad=("heterodyad_f05", "median"), n=("measurement_id", "size"))
    spread = med.groupby("canonical_pair_id").dyad.agg(
        n_sources="count", q05=lambda x: x.quantile(.05),
        median="median", q95=lambda x: x.quantile(.95)).reset_index()
    spread["width90"] = spread.q95 - spread.q05
    spread = spread[spread.n_sources > 1].copy()
    # Separate observed composition variation (pair-median RR, varying f) from
    # RR/source variation (citation medians, fixed f=0.5).
    comp_rows = []
    for pid, tgp in tg.groupby("canonical_pair_id"):
        rrp = rr[rr.canonical_pair_id == pid]
        fs = np.sort(tgp.fA.dropna().unique())
        if len(rrp) and len(fs) > 1:
            h = sequence_descriptors(np.repeat(rrp.r1.median(), len(fs)),
                                     np.repeat(rrp.r2.median(), len(fs)), fs)
            hv = np.asarray(h["dyad_12"]) + np.asarray(h["dyad_21"])
            comp_rows.append({"canonical_pair_id": pid,
                              "composition_width90": float(np.quantile(hv, .95)-np.quantile(hv, .05))})
    comp = pd.DataFrame(comp_rows)
    summary = {
        "n_rr_rows": int(len(rr)),
        "n_tg_pairs_with_rr": int(rr.canonical_pair_id.nunique()),
        "n_rr_rows_with_bibliographic_source": int((rr.citation_source != "unknown").sum()),
        "n_distinct_bibliographic_sources": int(rr.loc[rr.citation_source != "unknown", "citation_source"].nunique()),
        "n_source_pair_medians": int(len(med)),
        "n_pairs_with_multiple_sources": int((spread.n_sources > 1).sum()),
        "cross_citation_dyad_width90_at_f05_median": float(spread.width90.median()),
        "cross_citation_dyad_width90_at_f05_q75": float(spread.width90.quantile(.75)),
        "composition_only_dyad_width90_median": float(comp.composition_width90.median()),
        "composition_only_dyad_width90_q75": float(comp.composition_width90.quantile(.75)),
        "legacy_mixed_width90_median": 0.3198952942603245,
        "legacy_mixed_width90_q75": 0.4429223571312556,
        "legacy_width_note": "reproduced prior calculation pools source medians and each Tg sample composition in one distribution; it is not a source-only RR spread",
        "metadata_coverage": complete,
        "within_pair_condition_decomposition": decomp,
        "within_pair_and_citation_condition_decomposition": decomp_source_controlled,
        "interpretation": "each one-factor fraction is the reduction in within-stratum sum of squares explained by condition-cell means, using only replicated stratum-condition cells; factors overlap and shares are not additive",
        "limitations": [
            "condition fields are incomplete and strongly confounded with source",
            "source-controlled condition decomposition has low coverage because conditions rarely vary within a cited pair; inspect its row/group counts",
            "source aggregation can combine multiple conditions and fitting choices",
            "nominal temperature levels are rounded to 1 K; unreported conditions are an explicit category",
            "the decomposition is descriptive, not causal, and is not a measurement-error model",
        ],
    }
    spread.to_csv(PROJECT / "p4_rr_source_dyad_spread_f05.csv", index=False)
    return summary


def _metrics(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float)
    return {"n": int(len(y)), "MAE_K": float(np.mean(np.abs(y-p))),
            "RMSE_K": float(np.sqrt(np.mean((y-p)**2)))}


def _source_holdout(train: pd.DataFrame, test: pd.DataFrame) -> dict:
    models = {
        "rf_comp": (RandomForestRegressor(n_estimators=400, random_state=0, min_samples_leaf=2, n_jobs=-1), ["fA", "TgA_K", "TgB_K", "deltaTg_K"]),
        "rf_comp_seq": (RandomForestRegressor(n_estimators=400, random_state=0, min_samples_leaf=2, n_jobs=-1), ["fA", "TgA_K", "TgB_K", "deltaTg_K", "blockiness", "alternation", "seq_entropy", "dyad_11", "dyad_12", "dyad_21", "dyad_22"]),
    }
    ytr = train.Tg_exp_K.to_numpy() - train.Tg_Fox_w_K.to_numpy()
    out = {"fox": _metrics(test.Tg_exp_K, test.Tg_Fox_w_K)}
    for name, (model, cols) in models.items():
        Xtr = np.nan_to_num(train[cols].to_numpy(float), nan=0.0)
        Xte = np.nan_to_num(test[cols].to_numpy(float), nan=0.0)
        m = clone(model).fit(Xtr, ytr)
        pred = test.Tg_Fox_w_K.to_numpy() + m.predict(Xte)
        out[name] = _metrics(test.Tg_exp_K, pred)
        out[name]["delta_MAE_vs_Fox_K"] = out[name]["MAE_K"] - out["fox"]["MAE_K"]
    return out


def source_cv() -> dict:
    df = pd.read_parquet(PROJECT / "real_copolymer_tg_dataset.parquet")
    overlap = set(df[df.source == "tgvalue"].canonical_pair_id) & set(
        df[df.source == "styrenic"].canonical_pair_id)
    results = {}
    for held in sorted(df.source.unique()):
        train = df[df.source != held].copy()
        test = df[df.source == held].copy()
        overlap_here = set(test.canonical_pair_id) & set(train.canonical_pair_id)
        strict = test[~test.canonical_pair_id.isin(overlap_here)].copy()
        results[f"holdout_{held}"] = {
            "all_rows": _source_holdout(train, test),
            "strict_no_pair_overlap": _source_holdout(train, strict) if len(strict) else {},
            "n_train_rows": int(len(train)), "n_test_rows": int(len(test)),
            "n_test_pairs": int(test.canonical_pair_id.nunique()),
            "overlapping_pairs_removed_in_strict": int(len(overlap_here)),
            "strict_test_rows": int(len(strict)),
            "strict_test_pairs": int(strict.canonical_pair_id.nunique()),
        }
    results["overlap_pair_count"] = int(len(overlap))
    results["n_sources"] = int(df.source.nunique())
    results["interpretation"] = "only two Tg data sources are present; strict holdout removes any test monomer pair seen in training, so it is an intentionally small and noisy source-transfer check"
    return results


def main():
    out = {"rr_condition_audit": condition_spread(), "tg_source_holdout": source_cv()}
    (PROJECT / "p4_condition_source_cv.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
