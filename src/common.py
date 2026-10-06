"""Shared helpers for the reactivity-ratio research project.

Paths are resolved relative to the repository root so every script can be run
from anywhere with ``python -m`` or directly.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
STAGING = DATA / "staging"
STANDARDIZED = DATA / "standardized"
PROVENANCE = DATA / "provenance"
SPLITS = DATA / "benchmark_splits"
METADATA = ROOT / "metadata"
RESULTS = ROOT / "results"
PROJECTS = ROOT / "projects"

for _p in (STAGING, STANDARDIZED, PROVENANCE, SPLITS, METADATA, RESULTS,
           PROJECTS / "P1_rr_benchmark", PROJECTS / "P2_temperature_rr",
           PROJECTS / "P3_active_rr_design", PROJECTS / "P4_sequence_tg"):
    _p.mkdir(parents=True, exist_ok=True)

R_GAS = 8.31446261815324  # J / mol / K


def inchikey(smiles: str | float) -> str | None:
    """Return the first block (skeleton) of the InChIKey for a SMILES."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol).split("-")[0]
    except Exception:
        return None


def canonical_smiles(smiles: str | float) -> str | None:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def c_to_k(value):
    """Convert a Celsius value (or ``(min, nominal, max)``) to Kelvin scalars."""
    if value is None:
        return None
    return value + 273.15
