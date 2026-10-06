"""Separate RR-source spread from composition spread in predicted heterodyads.

The former S3 calculation pooled RR-source medians and the different sample
compositions in each pair into one distribution. This revision reports those
components separately. Cross-citation widths remain descriptive heterogeneity,
not formal uncertainty intervals.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import STANDARDIZED
from p4_sequence import sequence_descriptors

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projects" / "P4_sequence_tg"


def _width(values: np.ndarray) -> float:
    return float(np.quantile(values, .95) - np.quantile(values, .05))


def _citation_key(df: pd.DataFrame) -> pd.Series:
    ref = df.source_reference.fillna("").astype(str).str.strip()
    doi = df.source_doi.fillna("").astype(str).str.strip()
    valid_ref = ~ref.str.lower().isin(["", "na", "nan", "none"])
    valid_doi = ~doi.str.lower().isin(["", "na", "nan", "none"])
    key = ref.where(valid_ref, doi.where(valid_doi, ""))
    return key.str.lower().str.replace(r"\s+", " ", regex=True)


def main() -> None:
    tg = pd.read_parquet(PROJECT / "real_copolymer_tg_dataset.parquet")
    rr = pd.read_parquet(STANDARDIZED / "rrmaster_v1.parquet")
    rr = rr[rr.is_primary_measurement & (rr.r1 > 0) & (rr.r2 > 0)
            & rr.canonical_pair_id.isin(tg.canonical_pair_id.unique())].copy()
    rr["citation_key"] = _citation_key(rr)
    rr["heterodyad_f05"] = (
        np.asarray(sequence_descriptors(rr.r1, rr.r2, np.full(len(rr), .5))["dyad_12"])
        + np.asarray(sequence_descriptors(rr.r1, rr.r2, np.full(len(rr), .5))["dyad_21"]))

    rows = []
    for pair_id, tg_pair in tg.groupby("canonical_pair_id", sort=True):
        rp = rr[rr.canonical_pair_id == pair_id]
        if rp.empty:
            continue
        fs = np.sort(tg_pair.fA.dropna().unique())
        # Historical pooled width, retained to show what was confounded.
        ds = rp.groupby("source_dataset", as_index=False).agg(
            r1=("r1", "median"), r2=("r2", "median"))
        mix_desc = sequence_descriptors(np.repeat(ds.r1.to_numpy(), len(fs)),
                                        np.repeat(ds.r2.to_numpy(), len(fs)),
                                        np.tile(fs, len(ds)))
        mix = np.asarray(mix_desc["dyad_12"]) + np.asarray(mix_desc["dyad_21"])
        # Citation/source spread at common equimolar composition.
        cited = rp[rp.citation_key != ""].groupby("citation_key", as_index=False).agg(
            r1=("r1", "median"), r2=("r2", "median"))
        cite_desc = sequence_descriptors(cited.r1.to_numpy(), cited.r2.to_numpy(),
                                         np.full(len(cited), .5))
        cite_h = np.asarray(cite_desc["dyad_12"]) + np.asarray(cite_desc["dyad_21"])
        # Composition spread at fixed pair-median RR.
        comp_desc = sequence_descriptors(np.repeat(rp.r1.median(), len(fs)),
                                         np.repeat(rp.r2.median(), len(fs)), fs)
        comp_h = np.asarray(comp_desc["dyad_12"]) + np.asarray(comp_desc["dyad_21"])
        rows.append({
            "canonical_pair_id": pair_id,
            "nameA": tg_pair.nameA.iloc[0], "nameB": tg_pair.nameB.iloc[0],
            "n_rr_records": int(len(rp)), "n_citations": int(len(cited)),
            "n_sample_compositions": int(len(fs)),
            "legacy_mixed_width90": _width(mix),
            "citation_RR_width90_at_f05": _width(cite_h) if len(cited) > 1 else np.nan,
            "composition_only_width90_at_pair_median_RR": _width(comp_h) if len(fs) > 1 else 0.,
        })
    out_df = pd.DataFrame(rows)
    out_df.to_csv(PROJECT / "p4_supplementary_rr_spread.csv", index=False)
    summary = {
        "n_pairs": int(len(out_df)),
        "pairs_with_multiple_cited_sources": int(out_df.citation_RR_width90_at_f05.notna().sum()),
        "citation_RR_width90_at_f05_median": float(out_df.citation_RR_width90_at_f05.median()),
        "citation_RR_width90_at_f05_q75": float(out_df.citation_RR_width90_at_f05.quantile(.75)),
        "composition_only_width90_median": float(out_df.composition_only_width90_at_pair_median_RR.median()),
        "composition_only_width90_q75": float(out_df.composition_only_width90_at_pair_median_RR.quantile(.75)),
        "legacy_mixed_width90_median": float(out_df.legacy_mixed_width90.median()),
        "legacy_mixed_width90_q75": float(out_df.legacy_mixed_width90.quantile(.75)),
        "caution": "The legacy mixed width combines RR/source and sample composition variation. The two separated widths are not additive, and citation-source spread is not a confidence interval.",
    }
    (PROJECT / "p4_supplementary_rr_spread.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
