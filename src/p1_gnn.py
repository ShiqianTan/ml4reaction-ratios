"""P1 ReactivityRatioBench: GNN baseline (GIN on molecular graphs).

Predicts (ln r1, ln r2) jointly from the molecular graphs of the two monomers
using a shared GIN encoder + MLP head. Reuses the exact split schemes and
metric conventions of p1_benchmark.py so results are directly comparable.

Run: python3 src/p1_gnn.py [--quick]
Outputs: results/p1_gnn_results.csv, results/p1_gnn_oof.parquet
"""
from __future__ import annotations

import sys
import time
import numpy as np
import pandas as pd

from common import SPLITS, METADATA, RESULTS
from p1_benchmark import (random_split, group_split, one_monomer_out,
                          both_monomers_out, metrics, mayo_lewis_curve_error)

QUICK = "--quick" in sys.argv
HIDDEN = 128
EPOCHS = 25 if QUICK else 150
SEED = 0


# --------------------------------------------------------------------------- #
# graph construction (RDKit -> torch_geometric)
# --------------------------------------------------------------------------- #

ATOM_SYMS = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "Si", "B"]


def featurize_atom(a):
    from rdkit import Chem
    f = [1.0 if a.GetSymbol() == s else 0.0 for s in ATOM_SYMS]
    d = a.GetDegree()
    f += [1.0 if d == k else 0.0 for k in range(6)]
    f += [float(a.GetFormalCharge()), float(a.GetTotalNumHs()),
          float(a.GetIsAromatic()), float(a.IsInRing())]
    h = a.GetHybridization()
    f += [1.0 if h == t else 0.0 for t in
          [Chem.rdchem.HybridizationType.SP,
           Chem.rdchem.HybridizationType.SP2,
           Chem.rdchem.HybridizationType.SP3]]
    return f


def mol_to_graph(smiles):
    from rdkit import Chem
    import torch
    from torch_geometric.data import Data
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    xs = [featurize_atom(a) for a in mol.GetAtoms()]
    x = torch.tensor(xs, dtype=torch.float)
    src, dst = [], []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        src += [i, j]
        dst += [j, i]
    edge_index = torch.tensor([src, dst], dtype=torch.long) if src else \
        torch.zeros((2, 0), dtype=torch.long)
    return Data(x=x, edge_index=edge_index)


def build_graph_cache(reg):
    import torch
    cache = {}
    for ik, s in zip(reg.inchikey, reg.smiles_canonical):
        if ik not in cache and isinstance(s, str):
            g = mol_to_graph(s)
            if g is not None:
                cache[ik] = g
    return cache


# --------------------------------------------------------------------------- #
# model
# --------------------------------------------------------------------------- #

def build_model(node_dim):
    import torch
    import torch.nn as nn
    from torch_geometric.nn import GINConv, global_add_pool, global_mean_pool

    class GIN(nn.Module):
        def __init__(self):
            super().__init__()
            self.conv1 = GINConv(nn.Sequential(nn.Linear(node_dim, HIDDEN),
                                               nn.ReLU(), nn.Linear(HIDDEN, HIDDEN)))
            self.conv2 = GINConv(nn.Sequential(nn.Linear(HIDDEN, HIDDEN),
                                               nn.ReLU(), nn.Linear(HIDDEN, HIDDEN)))
            self.conv3 = GINConv(nn.Sequential(nn.Linear(HIDDEN, HIDDEN),
                                               nn.ReLU(), nn.Linear(HIDDEN, HIDDEN)))
            self.head = nn.Sequential(
                nn.Linear(2 * HIDDEN, HIDDEN), nn.ReLU(),
                nn.Linear(HIDDEN, 2))

        def encode(self, data):
            x, ei = data.x, data.edge_index
            x = self.conv1(x, ei).relu()
            x = self.conv2(x, ei).relu()
            x = self.conv3(x, ei).relu()
            return global_add_pool(x, data.batch) + global_mean_pool(x, data.batch)

        def forward(self, dataA, dataB):
            eA, eB = self.encode(dataA), self.encode(dataB)
            return self.head(torch.cat([eA, eB], dim=-1))

    return GIN()


def train_predict(model, cache, m1_ik, m2_ik, y, train_idx, test_idx):
    import torch
    from torch_geometric.data import Batch
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    lossf = torch.nn.MSELoss(reduction="none")

    def make_batch(rows):
        return (Batch.from_data_list([cache[m1_ik[i]] for i in rows]),
                Batch.from_data_list([cache[m2_ik[i]] for i in rows]))

    yt = torch.tensor(y[train_idx], dtype=torch.float)
    for _ in range(EPOCHS):
        model.train()
        opt.zero_grad()
        dA, dB = make_batch(train_idx)
        pred = model(dA, dB)
        loss = (lossf(pred, yt) * ~torch.isnan(yt)).sum() / \
            (~torch.isnan(yt)).sum()
        loss.backward()
        opt.step()

    model.eval()
    with torch.no_grad():
        dA, dB = make_batch(test_idx)
        return model(dA, dB).numpy()


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

def main():
    t0 = time.time()
    reg = pd.read_csv(METADATA / "monomer_registry.csv").drop_duplicates("inchikey")
    cache = build_graph_cache(reg)
    node_dim = len(ATOM_SYMS) + 6 + 4 + 3
    print(f"graph cache: {len(cache)} monomers, node_dim={node_dim}")

    df = pd.read_parquet(SPLITS / "dataset_C_condition.parquet").reset_index(drop=True)
    df = df[df.ln_r1.notna() & df.ln_r2.notna()].reset_index(drop=True)
    # only rows whose monomers are in the graph cache
    ok = df.m1_inchikey.isin(cache) & df.m2_inchikey.isin(cache)
    df = df[ok].reset_index(drop=True)
    print(f"benchmark rows: {len(df):,}")

    m1_ik = df.m1_inchikey.to_numpy()
    m2_ik = df.m2_inchikey.to_numpy()
    y = np.column_stack([df.ln_r1.to_numpy(float), df.ln_r2.to_numpy(float)])

    split_defs = {
        "random": random_split(df),
        "pair": group_split(df, "canonical_pair_id"),
        "one_monomer_out": one_monomer_out(df),
        "both_monomers_out": both_monomers_out(df),
    }

    import torch
    torch.manual_seed(SEED)
    np.random.seed(SEED)

    all_rows = []
    oof_store = {}
    for sname, splits in split_defs.items():
        if not splits:
            continue
        oof = np.full((len(df), 2), np.nan)
        for tr, te in splits:
            model = build_model(node_dim)
            oof[te] = train_predict(model, cache, m1_ik, m2_ik, y, tr, te)
            print(f"  {sname} fold done ({time.time()-t0:.0f}s)")
        oof_store[sname] = oof
        for j, target in enumerate(["ln_r1", "ln_r2"]):
            m = ~np.isnan(oof[:, j])
            all_rows.append({"target": target, "split": sname,
                             "features": "gnn", "model": "gin",
                             "n_eval": int(m.sum()),
                             **metrics(y[m, j], oof[m, j])})

    res = pd.DataFrame(all_rows)
    # E_curve from joint OOF
    r1t, r2t = df.r1.to_numpy(float), df.r2.to_numpy(float)
    for sname, oof in oof_store.items():
        m = ~np.isnan(oof[:, 0]) & ~np.isnan(oof[:, 1])
        e = mayo_lewis_curve_error(r1t[m], r2t[m], np.exp(oof[m, 0]),
                                   np.exp(oof[m, 1]))
        res.loc[(res.split == sname), "E_curve"] = e

    res.to_csv(RESULTS / "p1_gnn_results.csv", index=False)
    oof_df = df[["canonical_pair_id", "m1_inchikey", "m2_inchikey", "r1", "r2"]].copy()
    for sname, oof in oof_store.items():
        oof_df[f"{sname}_ln_r1"], oof_df[f"{sname}_ln_r2"] = oof[:, 0], oof[:, 1]
    oof_df.to_parquet(RESULTS / "p1_gnn_oof.parquet", index=False)

    print("\n=== GNN results ===")
    print(res[["target", "split", "MAE", "R2", "NLL", "coverage95", "E_curve"]]
          .round(4).to_string(index=False))
    print(f"\nwrote results ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
