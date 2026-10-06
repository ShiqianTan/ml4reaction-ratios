"""P4b: build homopolymer + alternating-copolymer Tg datasets.

Source: polyinfo_homopolymer_data.xlsx (PolyInfo) — 18k Tg entries including
~6.9k alternating (`-alt-`) copolymers.  Fox equation provides the composition
baseline; the residual is the sequence/interaction term.

Run: python3 src/p4_tg_data.py
"""
from __future__ import annotations

import json
import re
import numpy as np
import pandas as pd

from common import RAW, PROJECTS

P4 = PROJECTS / "P4_sequence_tg"
TG_COL = "Glass transition temperature_value_median"


def _strip_outer(s: str) -> str:
    s = s.strip()
    pairs = {"(": ")", "[": "]", "{": "}"}
    if s and s[0] in pairs and s[-1] == pairs[s[0]]:
        return s[1:-1]
    return s


def inner_name(name: str) -> str:
    s = str(name).strip()
    if s.lower().startswith("poly"):
        s = s[4:]
    return _strip_outer(s).strip()


def split_alt(inner: str):
    """Bracket-aware split on top-level '-alt-'."""
    depth = 0
    i = 0
    while i < len(inner):
        ch = inner[i]
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif depth == 0 and inner.startswith("-alt-", i):
            return _strip_outer(inner[:i]), _strip_outer(inner[i + 5:])
        i += 1
    return None


def load_polyinfo():
    x = pd.read_excel(RAW / "local/polyinfo_homopolymer_data.xlsx")
    x = x[["polymer_name", "SMILES", TG_COL]].dropna(subset=[TG_COL])
    x = x.rename(columns={TG_COL: "Tg_C", "SMILES": "smiles"})
    x["Tg_K"] = x.Tg_C + 273.15
    x["inner"] = x.polymer_name.map(inner_name)
    x = x.dropna(subset=["inner"]).drop_duplicates("inner")
    return x


def main():
    x = load_polyinfo()
    is_alt = x.polymer_name.astype(str).str.contains("-alt-", regex=False)
    mono = x[~is_alt & ~x.polymer_name.astype(str).str.contains("-co-", regex=False)]
    alt = x[is_alt].copy()
    print(f"polyinfo: {len(x)} unique-named Tg rows | "
          f"homo {len(mono)} | alt-copol {len(alt)}")

    # homopolymer lookup by inner name
    homo_lookup = mono.set_index("inner")[["Tg_C", "Tg_K", "smiles"]]
    mono_names = set(mono.inner)

    parts = alt.polymer_name.map(lambda n: split_alt(inner_name(n)))
    alt["comp1"] = [p[0] if p else None for p in parts]
    alt["comp2"] = [p[1] if p else None for p in parts]
    alt = alt.dropna(subset=["comp1", "comp2"])

    # match components to homopolymers
    def match(c):
        return c if c in mono_names else None

    alt["c1"] = alt.comp1.map(match)
    alt["c2"] = alt.comp2.map(match)
    matched = alt.dropna(subset=["c1", "c2"]).copy()
    print(f"alt copolymers with both components matched to a homopolymer: "
          f"{len(matched)} / {len(alt)}")

    matched["Tg1_C"] = homo_lookup.loc[matched.c1, "Tg_C"].to_numpy()
    matched["Tg2_C"] = homo_lookup.loc[matched.c2, "Tg_C"].to_numpy()
    matched["smiles1"] = homo_lookup.loc[matched.c1, "smiles"].to_numpy()
    matched["smiles2"] = homo_lookup.loc[matched.c2, "smiles"].to_numpy()

    # Fox baseline for strictly alternating (1:1 molar) copolymer, Kelvin
    T1 = matched.Tg1_C + 273.15
    T2 = matched.Tg2_C + 273.15
    matched["Fox_Tg_K"] = 1.0 / (0.5 / T1 + 0.5 / T2)
    matched["Fox_Tg_C"] = matched.Fox_Tg_K - 273.15
    matched["residual_C"] = matched.Tg_C - matched.Fox_Tg_C

    mono.to_parquet(P4 / "homopolymer_tg.parquet", index=False)
    matched.to_parquet(P4 / "alt_copolymer_tg.parquet", index=False)

    stats = {
        "n_homopolymer": int(len(mono)),
        "n_alt_copolymer_total": int(len(alt)),
        "n_alt_matched_components": int(len(matched)),
        "unique_component_pairs": int(
            matched[["c1", "c2"]].drop_duplicates().shape[0]),
        "unique_components": int(pd.unique(
            pd.concat([matched.c1, matched.c2])).size),
        "Tg_range_homo_C": [round(float(mono.Tg_C.min()), 1),
                            round(float(mono.Tg_C.max()), 1)],
        "Tg_range_copol_C": [round(float(matched.Tg_C.min()), 1),
                             round(float(matched.Tg_C.max()), 1)],
        "fox_MAE_C_vs_experiment": round(
            float(np.mean(np.abs(matched.residual_C))), 2),
        "fox_RMSE_C": round(float(np.sqrt(np.mean(matched.residual_C ** 2))), 2),
        "residual_std_C": round(float(matched.residual_C.std()), 2),
    }
    print(json.dumps(stats, indent=2))
    (P4 / "p4_data_summary.json").write_text(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
