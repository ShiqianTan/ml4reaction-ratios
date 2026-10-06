"""P1 Stage A-G: build the provenance-aware RRMaster v1 database.

Run:  python3 src/p1_build_rrmaster.py

Outputs (under data/ and metadata/):
  data/staging/all_raw_measurements.parquet
  data/standardized/rrmaster_v1.parquet
  metadata/monomer_registry.csv
  metadata/solvent_registry.csv
  metadata/source_lineage.csv
  data/provenance/duplicate_groups.csv
"""
from __future__ import annotations

import re
import numpy as np
import pandas as pd

from common import (RAW, STAGING, STANDARDIZED, PROVENANCE, METADATA,
                    inchikey, canonical_smiles, R_GAS)

R_EPS = 1e-9

# --------------------------------------------------------------------------- #
# Stage A: raw ingestion
# --------------------------------------------------------------------------- #

RAW_COLUMNS = [
    "source_dataset", "source_file", "source_row_id", "source_doi",
    "source_reference", "extraction_type",
    "monomer_1_name", "monomer_1_smiles", "monomer_2_name", "monomer_2_smiles",
    "r1_raw", "r2_raw", "r1_unc", "r2_unc",
    "temperature_raw", "temperature_unit",
    "solvent", "solvent_smiles", "initiator", "polymerization_method",
    "fitting_method", "notes",
]


def _blank(df: pd.DataFrame) -> pd.DataFrame:
    for c in RAW_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    return df[RAW_COLUMNS]


def load_chen_jackson_handbook() -> pd.DataFrame:
    smi = RAW / "chen_jackson/dataset/monomer_react_ratio_smarts.csv"
    err = RAW / "chen_jackson/dataset/Handbook_original_filtered.csv"
    df = pd.read_csv(smi)
    df["r1"] = pd.to_numeric(df["r1"], errors="coerce")
    df["r2"] = pd.to_numeric(df["r2"], errors="coerce")
    df = df.rename(columns={"Monomer_1": "monomer_1_name",
                            "Monomer_1_smi": "monomer_1_smiles",
                            "Monomer_2": "monomer_2_name",
                            "Monomer_2_smi": "monomer_2_smiles",
                            "r1": "r1_raw", "r2": "r2_raw"})
    # attach uncertainty from the aligned handbook table (match on names+rr)
    e = pd.read_csv(err)
    e["key"] = (e.Monomer_1.astype(str) + "|" + e.Monomer_2.astype(str) + "|"
                + e.r1.round(4).astype(str) + "|" + e.r2.round(4).astype(str))
    e = e.drop_duplicates("key").set_index("key")
    df["key"] = (df.monomer_1_name.astype(str) + "|" + df.monomer_2_name.astype(str)
                 + "|" + df.r1_raw.round(4).astype(str) + "|" + df.r2_raw.round(4).astype(str))
    df["r1_unc"] = df.key.map(e.r1_err)
    df["r2_unc"] = df.key.map(e.r2_err)
    df["source_dataset"] = "chen_jackson_handbook"
    df["source_file"] = smi.name
    df["source_row_id"] = np.arange(len(df))
    df["extraction_type"] = "OCR/handbook"
    df["notes"] = "Polymer Handbook RR table (secondary source)"
    return _blank(df)


def _load_chen_jackson_aux(fname: str, dataset: str) -> pd.DataFrame:
    p = RAW / "chen_jackson/dataset" / fname
    df = pd.read_csv(p)
    df = df.rename(columns={"Monomer_1": "monomer_1_name",
                            "Monomer_1_smi": "monomer_1_smiles",
                            "Monomer_2": "monomer_2_name",
                            "Monomer_2_smi": "monomer_2_smiles",
                            "r1": "r1_raw", "r2": "r2_raw",
                            "Temperature (K)": "temperature_raw"})
    df["temperature_unit"] = "K"
    df["source_dataset"] = dataset
    df["source_file"] = fname
    df["source_row_id"] = np.arange(len(df))
    df["extraction_type"] = "original_csv"
    return _blank(df)


def load_polycarp() -> pd.DataFrame:
    p = RAW / "polycarp/copol_prediction/processed_data.csv"
    df = pd.read_csv(p)
    out = pd.DataFrame({
        "source_dataset": "polycarp",
        "source_file": p.name,
        "source_row_id": np.arange(len(df)),
        "source_doi": df.get("original_source"),
        "source_reference": df.get("PDF_name"),
        "extraction_type": "VLM",
        "monomer_1_name": df.get("monomer1_name"),
        "monomer_1_smiles": df.get("monomer1_smiles"),
        "monomer_2_name": df.get("monomer2_name"),
        "monomer_2_smiles": df.get("monomer2_smiles"),
        "r1_raw": pd.to_numeric(df.get("constant_1"), errors="coerce"),
        "r2_raw": pd.to_numeric(df.get("constant_2"), errors="coerce"),
        "r1_unc": pd.to_numeric(df.get("constant_conf_1"), errors="coerce"),
        "r2_unc": pd.to_numeric(df.get("constant_conf_2"), errors="coerce"),
        "temperature_raw": pd.to_numeric(df.get("temperature"), errors="coerce"),
        "temperature_unit": df.get("temperature_unit"),
        "solvent": df.get("solvent"),
        "solvent_smiles": df.get("solvent_smiles"),
        "initiator": np.nan,
        "polymerization_method": df.get("method"),
        "fitting_method": df.get("determination_method"),
        "notes": df.get("polymerization_type"),
    })
    return _blank(out)


def load_coxwell_full() -> pd.DataFrame:
    p = RAW / "coxwell/Coxwell_RR full dataset final.xlsx"
    df = pd.read_excel(p)
    out = pd.DataFrame({
        "source_dataset": "coxwell_full",
        "source_file": p.name,
        "source_row_id": np.arange(len(df)),
        "source_doi": np.nan,
        "source_reference": df.get("Source"),
        "extraction_type": "SI_excel",
        "monomer_1_name": df.get("monomer 1"),
        "monomer_1_smiles": df.get("M1 Smiles"),
        "monomer_2_name": df.get("monomer 2"),
        "monomer_2_smiles": df.get("M2 Smiles"),
        "r1_raw": pd.to_numeric(df.get("r1"), errors="coerce"),
        "r2_raw": np.nan,  # only forward direction tabulated
        "r1_unc": np.nan,
        "r2_unc": np.nan,
        "temperature_raw": pd.to_numeric(df.get("Temp (K)"), errors="coerce"),
        "temperature_unit": "K",
        "solvent": df.get("Solvent"),
        "solvent_smiles": np.nan,
        "initiator": np.nan,
        "polymerization_method": np.nan,
        "fitting_method": df.get("Stat. Method"),
        "notes": df.get("Notes "),
    })
    # recover the reverse-direction r2 from symmetric pair rows
    out = _recover_symmetric_r2(out)
    return _blank(out)


def load_coxwell_error() -> pd.DataFrame:
    p = RAW / "coxwell/Coxwell_RR_error.xlsx"
    df = pd.read_excel(p, sheet_name="Data")
    out = pd.DataFrame({
        "source_dataset": "coxwell_error",
        "source_file": p.name,
        "source_row_id": np.arange(len(df)),
        "source_doi": np.nan,
        "source_reference": df.get("Source"),
        "extraction_type": "SI_excel",
        "monomer_1_name": df.get("M1"),
        "monomer_1_smiles": np.nan,
        "monomer_2_name": df.get("M2"),
        "monomer_2_smiles": np.nan,
        "r1_raw": pd.to_numeric(df.get("r1"), errors="coerce"),
        "r2_raw": pd.to_numeric(df.get("r2"), errors="coerce"),
        "r1_unc": np.nan,
        "r2_unc": np.nan,
        "temperature_raw": pd.to_numeric(df.get("Temperature (K)"), errors="coerce"),
        "temperature_unit": "K",
        "solvent": df.get("Solvent"),
        "solvent_smiles": np.nan,
        "initiator": np.nan,
        "polymerization_method": np.nan,
        "fitting_method": df.get("Method"),
        "notes": "multi-temperature / uncertainty study",
    })
    return _blank(out)


def _recover_symmetric_r2(df: pd.DataFrame) -> pd.DataFrame:
    """Fill r2 for a row A-B using the r1 of the tabulated B-A row."""
    key = df["monomer_1_name"].astype(str) + "||" + df["monomer_2_name"].astype(str)
    rev = df.assign(_revkey=df["monomer_2_name"].astype(str) + "||"
                    + df["monomer_1_name"].astype(str))
    lookup = dict(zip(rev["_revkey"], rev["r1_raw"]))
    df["r2_raw"] = df["r2_raw"].where(df["r2_raw"].notna(),
                                      key.map(lookup))
    return df


def stage_a() -> pd.DataFrame:
    frames = [
        load_chen_jackson_handbook(),
        _load_chen_jackson_aux("CKAs.csv", "chen_jackson_cka"),
        _load_chen_jackson_aux("textbook_monomer_react_ratio_smarts.csv",
                               "chen_jackson_textbook"),
        load_polycarp(),
        load_coxwell_full(),
        load_coxwell_error(),
    ]
    allraw = pd.concat(frames, ignore_index=True)
    allraw["retrieval_date"] = "2026-10-05"
    allraw.to_parquet(STAGING / "all_raw_measurements.parquet", index=False)
    print(f"[Stage A] raw measurements: {len(allraw):,}")
    print(allraw.source_dataset.value_counts().to_string())
    return allraw


# --------------------------------------------------------------------------- #
# Stage B: monomer identity standardisation
# --------------------------------------------------------------------------- #

def stage_b(df: pd.DataFrame):
    registry = {}

    def register(row_idx, name, smiles):
        ik = inchikey(smiles)
        csm = canonical_smiles(smiles)
        registry[(str(name), str(smiles))] = (ik, csm)
        return ik, csm

    m1 = [register(i, n, s) for i, (n, s) in
          enumerate(zip(df.monomer_1_name, df.monomer_1_smiles))]
    m2 = [register(i, n, s) for i, (n, s) in
          enumerate(zip(df.monomer_2_name, df.monomer_2_smiles))]
    df["m1_inchikey"] = [x[0] for x in m1]
    df["m1_canonical_smiles"] = [x[1] for x in m1]
    df["m2_inchikey"] = [x[0] for x in m2]
    df["m2_canonical_smiles"] = [x[1] for x in m2]

    df["smiles_valid"] = df.m1_inchikey.notna() & df.m2_inchikey.notna()
    # monomer registry
    reg_rows = []
    for key, (ik, csm) in registry.items():
        name, smi = key
        reg_rows.append({"monomer_name_raw": name, "smiles_raw": smi,
                         "smiles_canonical": csm, "inchikey": ik})
    reg = pd.DataFrame(reg_rows).drop_duplicates(subset=["smiles_raw"])
    reg["monomer_id"] = ["M%04d" % i for i in range(len(reg))]
    reg["manual_review_flag"] = reg.inchikey.isna()
    reg.to_csv(METADATA / "monomer_registry.csv", index=False)
    print(f"[Stage B] monomer registry: {len(reg):,} unique raw SMILES; "
          f"invalid={int(reg.inchikey.isna().sum())}")
    return df, reg


# --------------------------------------------------------------------------- #
# Stage C: RR direction canonicalisation
# --------------------------------------------------------------------------- #

def stage_c(df: pd.DataFrame) -> pd.DataFrame:
    order = df.m1_inchikey < df.m2_inchikey
    swap = (~order) & df.m1_inchikey.notna() & df.m2_inchikey.notna()
    for a, b in [("monomer_1_name", "monomer_2_name"),
                 ("monomer_1_smiles", "monomer_2_smiles"),
                 ("m1_inchikey", "m2_inchikey"),
                 ("m1_canonical_smiles", "m2_canonical_smiles"),
                 ("r1_raw", "r2_raw"), ("r1_unc", "r2_unc")]:
        df[a], df[b] = df[a].where(~swap, df[b]), df[b].where(~swap, df[a])
    df["orientation_swapped"] = swap
    df["canonical_pair_id"] = df.m1_inchikey + "__" + df.m2_inchikey
    df.loc[df.m1_inchikey.isna() | df.m2_inchikey.isna(), "canonical_pair_id"] = np.nan
    print(f"[Stage C] orientation swapped: {int(swap.sum()):,}")
    return df


# --------------------------------------------------------------------------- #
# Stage D: reaction-condition standardisation
# --------------------------------------------------------------------------- #

_SOLVENT_FAMILY = {
    "water": "water", "methanol": "alcohol", "ethanol": "alcohol",
    "benzene": "aromatic", "toluene": "aromatic", "xylene": "aromatic",
    "chloroform": "halogenated", "dichloromethane": "halogenated",
    "dimethylformamide": "amide", "n,n-dimethylformamide": "amide",
    "dmso": "sulfoxide", "dimethyl sulfoxide": "sulfoxide",
    "tetrahydrofuran": "ether", "thf": "ether", "dioxane": "ether",
    "acetonitrile": "nitrile", "bulk": "bulk", "none": "bulk",
    "hexane": "alkane", "cyclohexane": "alkane", "acetone": "ketone",
}


def _solvent_family(name):
    if not isinstance(name, str):
        return None
    s = name.strip().lower()
    if "bulk" in s:
        return "bulk"
    for k, v in _SOLVENT_FAMILY.items():
        if k in s:
            return v
    return "other"


def stage_d(df: pd.DataFrame) -> pd.DataFrame:
    t = pd.to_numeric(df.temperature_raw, errors="coerce")
    unit = df.temperature_unit.astype(str).str.lower()
    kelvin = t.where(unit.str.contains("k") & ~unit.str.contains("°c"), t + 273.15)
    df["temperature_K"] = np.where(t.notna(), kelvin, np.nan)
    df["temperature_K"] = pd.to_numeric(df.temperature_K, errors="coerce")
    # physical sanity window
    bad = (df.temperature_K < 150) | (df.temperature_K > 600)
    df.loc[bad, "temperature_K"] = np.nan
    df["temperature_nominal_K"] = df.temperature_K
    df["solvent"] = df.solvent.replace({np.nan: np.nan})
    df["solvent_family"] = df.solvent.map(_solvent_family)
    df["solvent_family"] = df.solvent_family.fillna("unknown")
    # solvent registry
    sreg = (df[["solvent", "solvent_smiles"]].dropna(subset=["solvent"])
            .drop_duplicates(subset=["solvent"]))
    sreg = sreg.assign(solvent_id=["S%04d" % i for i in range(len(sreg))],
                       family=sreg.solvent.map(_solvent_family))
    sreg.to_csv(METADATA / "solvent_registry.csv", index=False)
    print(f"[Stage D] temperature known: {int(df.temperature_K.notna().sum()):,}; "
          f"solvent known: {int((df.solvent_family != 'unknown').sum()):,}")
    return df


# --------------------------------------------------------------------------- #
# Stage E: RR targets / QC
# --------------------------------------------------------------------------- #

def stage_e(df: pd.DataFrame) -> pd.DataFrame:
    for c in ("r1_raw", "r2_raw"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["r1"] = df.r1_raw
    df["r2"] = df.r2_raw
    df["ln_r1"] = np.log(df.r1.where(df.r1 > R_EPS))
    df["ln_r2"] = np.log(df.r2.where(df.r2 > R_EPS))
    df["log10_r1"] = np.log10(df.r1.where(df.r1 > R_EPS))
    df["log10_r2"] = np.log10(df.r2.where(df.r2 > R_EPS))
    df["nonpositive_r1"] = df.r1 <= 0
    df["nonpositive_r2"] = df.r2 <= 0
    df["extreme_rr_flag"] = ((df.r1 > 100) | (df.r1 < 0.01)
                             | (df.r2 > 100) | (df.r2 < 0.01))
    rt = R_GAS * df.temperature_K / 1000.0  # kJ/mol
    df["ddg1_kJ_mol"] = -rt * df.ln_r1
    df["ddg2_kJ_mol"] = -rt * df.ln_r2
    both = df.r1.notna() & df.r2.notna()
    print(f"[Stage E] usable RR pairs (both numeric): {int(both.sum()):,}; "
          f"nonpositive r1={int(df.nonpositive_r1.sum())}, "
          f"r2={int(df.nonpositive_r2.sum())}")
    return df


# --------------------------------------------------------------------------- #
# Stage F: provenance & de-duplication
# --------------------------------------------------------------------------- #

_SOURCE_PRIORITY = {
    "coxwell_error": 1, "coxwell_full": 2, "polycarp": 3,
    "chen_jackson_cka": 4, "chen_jackson_textbook": 5,
    "chen_jackson_handbook": 6,
}


def stage_f(df: pd.DataFrame) -> pd.DataFrame:
    df = df.reset_index(drop=True)
    df["measurement_id"] = df.source_dataset + ":" + df.source_row_id.astype(str)
    df["_r1r"] = df.r1_raw.round(2)
    df["_r2r"] = df.r2_raw.round(2)
    df["_tt"] = df.temperature_nominal_K.round(-1)  # nearest 10 K
    df["_loose_key"] = (df.canonical_pair_id + "|"
                        + df._r1r.astype(str) + "|" + df._r2r.astype(str))
    df["_exact_key"] = (df.canonical_pair_id + "|"
                        + df.r1_raw.round(4).astype(str) + "|"
                        + df.r2_raw.round(4).astype(str))
    df["source_priority"] = df.source_dataset.map(_SOURCE_PRIORITY).fillna(99)

    # duplicate groups: rows sharing pair + rounded rr (any source)
    valid = df.canonical_pair_id.notna() & df.r1_raw.notna() & df.r2_raw.notna()
    codes, uniques = pd.factorize(df.loc[valid, "_loose_key"])
    grp = pd.Series(pd.NA, index=df.index, dtype="object")
    grp.loc[valid] = ["D%06d" % (c + 1) for c in codes]
    df["duplicate_group_id"] = grp

    # choose the most primary member per duplicate group
    vg = df[df.duplicate_group_id.notna()]
    idx = vg.groupby("duplicate_group_id")["source_priority"].idxmin()
    prim = vg.loc[idx, ["duplicate_group_id", "measurement_id", "source_dataset"]]
    map_id = dict(zip(prim.duplicate_group_id, prim.measurement_id))
    map_src = dict(zip(prim.duplicate_group_id, prim.source_dataset))
    df["primary_measurement_id"] = df.duplicate_group_id.map(map_id)
    df["primary_source"] = df.duplicate_group_id.map(map_src)
    df["is_primary_measurement"] = df.measurement_id == df.primary_measurement_id

    # secondary sources: all other sources present in the group
    def _others(s):
        return sorted(set(s) - {map_src.get(s.name)})
    sec = vg.groupby("duplicate_group_id")["source_dataset"].apply(_others)
    df["secondary_sources"] = df.duplicate_group_id.map(sec)
    nanmask = df.duplicate_group_id.isna()
    df.loc[nanmask, "is_primary_measurement"] = True
    df.loc[nanmask, "primary_measurement_id"] = df.loc[nanmask, "measurement_id"]
    df.loc[nanmask, "primary_source"] = df.loc[nanmask, "source_dataset"]
    df["secondary_sources"] = df.secondary_sources.apply(
        lambda x: x if isinstance(x, list) else [])

    dup = df[["measurement_id", "canonical_pair_id", "duplicate_group_id",
              "primary_measurement_id", "is_primary_measurement",
              "primary_source", "source_dataset"]].copy()
    dup.to_csv(PROVENANCE / "duplicate_groups.csv", index=False)

    lineage = (df.groupby("source_dataset")
               .agg(n=("measurement_id", "size"),
                    n_primary=("is_primary_measurement", "sum"),
                    n_doi=("source_doi", lambda x: x.notna().sum()))
               .reset_index())
    lineage.to_csv(METADATA / "source_lineage.csv", index=False)
    print(f"[Stage F] duplicate groups: {df.duplicate_group_id.nunique():,}; "
          f"primary measurements: {int(df.is_primary_measurement.sum()):,}")
    return df


# --------------------------------------------------------------------------- #
# Stage G: experimental quality score
# --------------------------------------------------------------------------- #

_MODERN_FIT = {"nlls", "nonlinear", "evm", "error-in-variables", "ymlr",
               "mayo-lewis nonlinear", "non-linear"}


def stage_g(df: pd.DataFrame) -> pd.DataFrame:
    traceable = df.source_doi.notna() | df.source_reference.notna()
    conditions = df.temperature_K.notna() & (df.solvent_family != "unknown")
    fm = df.fitting_method.fillna("").astype(str).str.lower()
    modern = fm.apply(lambda s: any(k in s for k in _MODERN_FIT))
    unc = df.r1_unc.notna() | df.r2_unc.notna()
    replicate = df.duplicate_group_id.duplicated(keep=False)
    score = (traceable.astype(int) * 20 + conditions.astype(int) * 15
             + modern.astype(int) * 15 + unc.astype(int) * 15
             + replicate.astype(int) * 10)
    # secondary handbook sources lose the traceability credit to primary papers
    score = score.where(df.source_dataset != "chen_jackson_handbook", score - 10)
    df["quality_score"] = score.clip(lower=0)
    df["quality_grade"] = pd.cut(df.quality_score, bins=[-1, 20, 40, 65, 100],
                                 labels=["D", "C", "B", "A"])
    print("[Stage G] quality grade distribution:")
    print(df.quality_grade.value_counts().sort_index().to_string())
    return df


def main():
    raw = stage_a()
    df, reg = stage_b(raw)
    df = stage_c(df)
    df = stage_d(df)
    df = stage_e(df)
    df = stage_f(df)
    df = stage_g(df)
    # drop scratch columns
    drop = [c for c in df.columns if c.startswith("_")]
    df = df.drop(columns=drop)
    out = STANDARDIZED / "rrmaster_v1.parquet"
    df.to_parquet(out, index=False)
    print(f"\n[RRMaster v1] rows={len(df):,}  -> {out}")
    print(f"  unique canonical pairs: {df.canonical_pair_id.nunique():,}")
    print(f"  unique monomers: "
          f"{pd.unique(pd.concat([df.m1_inchikey, df.m2_inchikey]).dropna()).size:,}")
    print(f"  usable both-direction RR: "
          f"{int((df.ln_r1.notna() & df.ln_r2.notna()).sum()):,}")


if __name__ == "__main__":
    main()
