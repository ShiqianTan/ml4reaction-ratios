"""P4c: real Tg models from repeat-unit structure.

  * homopolymer Tg      (11,380 polymer structures)
  * alternating-copolymer Tg (6,891 structures)

Splits: random vs scaffold-cluster (chemistry-aware OOD).  Models: Ridge / RF /
XGBoost.  Target: Tg (deg C).

Run: python3 src/p4_tg_model.py
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import PROJECTS

P4 = PROJECTS / "P4_sequence_tg"
NBITS = 2048


def morgan(smiles, nbits=NBITS):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    arr = np.zeros(nbits, np.int8)
    mol = Chem.MolFromSmiles(str(smiles)) if smiles is not None else None
    if mol is None:
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nbits)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def rdkit_block(smiles):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(str(smiles)) if smiles is not None else None
    names = ["MolWt", "MolLogP", "TPSA", "NumRotatableBonds", "FractionCSP3",
             "RingCount", "NumAromaticRings", "NumHDonors", "NumHAcceptors",
             "HeavyAtomCount"]
    out = []
    for n in names:
        try:
            f = getattr(rdMolDescriptors, n, None) or getattr(Descriptors, n)
            out.append(float(f(mol)) if mol is not None else np.nan)
        except Exception:
            out.append(np.nan)
    return np.array(out, np.float32)


def scaffold(smiles):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    mol = Chem.MolFromSmiles(str(smiles)) if smiles is not None else None
    if mol is None:
        return "NA"
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return "NA"


def features(smiles_list):
    X = np.vstack([np.concatenate([morgan(s), rdkit_block(s)])
                   for s in smiles_list])
    return np.nan_to_num(X, nan=0.0)


def evaluate(df, tag):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold, GroupKFold
    from sklearn.metrics import r2_score
    y = df.Tg_C.to_numpy(float)
    X = features(df.smiles)
    print(f"\n[{tag}] n={len(df)}  features={X.shape[1]}")
    out = {}
    folds = {"random": list(KFold(5, shuffle=True, random_state=0).split(X)),
             "scaffold": list(GroupKFold(5).split(
                 X, groups=df["_scaffold"].values))}
    for sname, sp in folds.items():
        for mname, model in [("ridge", Ridge(alpha=1.0)),
                             ("rf", RandomForestRegressor(
                                 n_estimators=300, min_samples_leaf=2,
                                 n_jobs=-1, random_state=0))]:
            from sklearn.base import clone
            oof = np.full(len(df), np.nan)
            for tr, te in sp:
                m = clone(model)
                m.fit(X[tr], y[tr])
                oof[te] = m.predict(X[te])
            met = {"MAE_C": round(float(np.mean(np.abs(y - oof))), 2),
                   "RMSE_C": round(float(np.sqrt(np.mean((y - oof) ** 2))), 2),
                   "R2": round(float(r2_score(y, oof)), 3)}
            out[f"{sname}_{mname}"] = met
            print(f"  {sname:9s} {mname:5s} MAE={met['MAE_C']:6.2f} C  "
                  f"R2={met['R2']:.3f}")
    return out


def main():
    from sklearn.model_selection import train_test_split
    homo = pd.read_parquet(P4 / "homopolymer_tg.parquet")
    alt = pd.read_parquet(P4 / "alt_copolymer_tg.parquet")

    # limit size for runtime; keep all alt, sample homo
    if len(homo) > 8000:
        homo = homo.sample(8000, random_state=0)
    for name, d in [("homo", homo), ("alt", alt)]:
        d["_scaffold"] = d.smiles.map(scaffold)

    res = {}
    res["homopolymer"] = evaluate(homo, "homopolymer")
    res["alt_copolymer"] = evaluate(alt, "alt-copolymer")

    (P4 / "p4_tg_model_results.json").write_text(json.dumps(res, indent=2))
    print("\nsaved", P4 / "p4_tg_model_results.json")


if __name__ == "__main__":
    main()
