"""Feature generation for the P1 ReactivityRatioBench.

  * ECFP4 (2048 bit) per monomer
  * RDKit 2D descriptors per monomer
  * condition features (temperature, solvent family, polymerisation method)
  * CopDDB physics descriptors (Delta-barrier, Delta-reaction) for r-direction
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import RAW, STANDARDIZED, METADATA, inchikey, canonical_smiles

ECFP_RADIUS = 2
ECFP_NBITS = 2048

_RDKIT_DESC = [
    "MolWt", "MolLogP", "MolMR", "TPSA", "NumHAcceptors", "NumHDonors",
    "NumRotatableBonds", "RingCount", "NumAromaticRings", "FractionCSP3",
    "HeavyAtomCount", "NumValenceElectrons", "MaxPartialCharge",
    "MinPartialCharge", "LabuteASA", "BalabanJ", "BertzCT", "HallKierAlpha",
    "Ipc", "Kappa1", "Kappa2", "Kappa3", "Chi0", "Chi1", "Chi2n", "Chi3n",
    "Chi4n", "Chi0v", "Chi1v", "Chi2v", "Chi3v", "Chi4v",
]


def _morgan(smiles, nbits=ECFP_NBITS, radius=ECFP_RADIUS):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    arr = np.zeros(nbits, dtype=np.int8)
    if mol is None:
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def _rdkit_descriptors(smiles):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    vals = []
    for name in _RDKIT_DESC:
        try:
            if mol is None:
                vals.append(np.nan)
            elif hasattr(rdMolDescriptors, name):
                vals.append(float(getattr(rdMolDescriptors, name)(mol)))
            else:
                vals.append(float(getattr(Descriptors, name)(mol)))
        except Exception:
            vals.append(np.nan)
    return np.array(vals, dtype=np.float32)


def monomer_features():
    """Return (ecfp_df, desc_df, smiles_by_inchikey) keyed by inchikey."""
    reg = pd.read_csv(METADATA / "monomer_registry.csv")
    reg = reg.dropna(subset=["inchikey"]).drop_duplicates("inchikey")
    ecfp, desc, smi = {}, {}, {}
    for ik, s in zip(reg.inchikey, reg.smiles_canonical):
        ecfp[ik] = _morgan(s)
        desc[ik] = _rdkit_descriptors(s)
        smi[ik] = s
    e = pd.DataFrame(ecfp).T
    e.columns = [f"ecfp_{i}" for i in range(e.shape[1])]
    d = pd.DataFrame(desc).T
    d.columns = _RDKIT_DESC
    return e, d, smi


# --------------------------------------------------------------------------- #
# CopDDB physics descriptors
# --------------------------------------------------------------------------- #

def _copddb_table():
    p = RAW / "copddb/copddb/datasets/data/PropagationQuantumChem_2024-03-12.csv"
    c = pd.read_csv(p)
    c["rad_ik"] = [inchikey(s) for s in c.Radical]
    c["mon_ik"] = [inchikey(s) for s in c.Monomer]
    return c


def copddb_rr_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    """Delta descriptors for r1 and r2 directions from CopDDB barriers.

    r1 = k11/k12 uses (R1,M1) vs (R1,M2):
        d_barrier_1 = B(R1,M2) - B(R1,M1)
    r2 = k22/k21 uses (R2,M2) vs (R2,M1):
        d_barrier_2 = B(R2,M1) - B(R2,M2)
    """
    c = _copddb_table()
    b = c.set_index(["rad_ik", "mon_ik"])[["DE_barrier", "DE_reaction",
                                           "DE_SHgap", "DE_SLgap"]]
    out = pd.DataFrame(index=df.index)
    keys_11 = pd.MultiIndex.from_arrays([df.m1_inchikey, df.m1_inchikey])
    keys_12 = pd.MultiIndex.from_arrays([df.m1_inchikey, df.m2_inchikey])
    keys_22 = pd.MultiIndex.from_arrays([df.m2_inchikey, df.m2_inchikey])
    keys_21 = pd.MultiIndex.from_arrays([df.m2_inchikey, df.m1_inchikey])

    def g(keys, col):
        return b[col].reindex(keys).to_numpy()

    for col in ["DE_barrier", "DE_reaction", "DE_SHgap", "DE_SLgap"]:
        b11, b12 = g(keys_11, col), g(keys_12, col)
        b22, b21 = g(keys_22, col), g(keys_21, col)
        out[f"copddb_d1_{col}"] = b12 - b11
        out[f"copddb_d2_{col}"] = b21 - b22
    out["copddb_available"] = out[[c for c in out if c.startswith("copddb_d1_")]].notna().any(axis=1)
    return out


def condition_features(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    out["temperature_K"] = df.temperature_K.astype(float)
    out["inv_temperature"] = 1000.0 / df.temperature_K
    fam = df.solvent_family.fillna("unknown")
    fams = ["bulk", "aromatic", "alcohol", "halogenated", "amide", "ether",
            "water", "alkane", "nitrile", "sulfoxide", "ketone", "other",
            "unknown"]
    for f in fams:
        out[f"solv_{f}"] = (fam == f).astype(int)
    return out
