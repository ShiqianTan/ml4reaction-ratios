# ML for Reaction Ratios — Reactivity-Ratio Research Program

Implementation of the four-project plan in
[`reactivity_ratio_project_plan.md`](reactivity_ratio_project_plan.md):

| Priority | Project | Status |
|---|---|---|
| **P1** | ReactivityRatioBench (data + benchmark) | 🟢 data base + benchmark v1 done |
| **P2** | Temperature-dependent RR thermodynamics | 🟢 done |
| **P3** | ML prior + Bayesian / active experimental design | 🟢 done |
| **P4** | RR → Sequence → Tg multi-scale Δ-learning | 🟢 done (real Tg data) |

> Progress log and key results are maintained below. Data lineage is frozen in
> `data/raw/` (never overwritten); all derived artefacts live in
> `data/staging`, `data/standardized`, `metadata`, `results`.

---

## Repository layout

```text
data/raw/            frozen source datasets (git clones + local copies)
data/staging/        all_raw_measurements.parquet (Stage A)
data/standardized/   rrmaster_v1.parquet (Stage G)
data/provenance/     duplicate_groups.csv
data/benchmark_splits/ dataset_U / HQ / C / legacy parquet
metadata/            monomer_registry.csv, solvent_registry.csv, source_lineage.csv
src/                 pipeline + benchmark scripts
results/             benchmark tables and figures
projects/            per-project outputs (P1..P4)
```

## Environment / how to run

Uses the base conda Python 3.11 (rdkit, pandas, scikit-learn, xgboost,
lightgbm, torch, pyarrow).

```bash
python3 src/p1_build_rrmaster.py     # Stage A-G  -> RRMaster v1
python3 src/p1_feasibility.py        # Phase-2 stats + dataset versions
python3 src/p1_benchmark.py          # P1 experiment matrix (use --quick to smoke-test)
```

---

## P1 — ReactivityRatioBench

### Data sources ingested (Stage A)

| source | rows | conditions | provenance |
|---|---:|---|---|
| PolyCARP (`processed_data.csv`) | 7,582 | T, solvent, method, Q-e, DOI | VLM extraction |
| Chen & Jackson / Polymer Handbook | 4,787 | — (r_err only) | handbook OCR (secondary) |
| Coxwell 2026 full dataset | 484 | T, solvent, stat. method, source | primary literature |
| Coxwell 2026 error study | 50 | T, solvent, method | multi-T / uncertainty |
| Chen & Jackson textbook | 79 | T | primary |
| Chen & Jackson CKAs | 10 | T | primary |
| **total raw** | **12,992** | | |

Physics descriptors: **CopDDB** (2,500 radical–monomer pairs, barrier/reaction
energies, orbital gaps, %Vbur). IUPAC `f0XFtestset.csv` reserved for P3.
Homopolymer + copolymer Tg data reserved for P4.

### Pipeline (Stages A–G)

* **B** monomer identity → canonical SMILES + InChIKey (**1,718** raw SMILES, 1 invalid).
* **C** RR direction canonicalisation (sorted InChIKey; 6,422 rows swapped).
* **D** conditions → `temperature_K` (K / °C unified, physical window 150–600 K),
  solvent registry + family (BULK kept as its own token).
* **E** targets `ln_r`, `log10_r`, `ddg_kJ_mol = -RT ln r`; QC flags for
  non-positive and extreme RR.
* **F** provenance de-duplication (pair + rounded r → duplicate groups; most
  primary source chosen; **302** groups span >1 source).
* **G** experimental quality score / grade A–D.

Output: **`data/standardized/rrmaster_v1.parquet`** — 12,992 rows,
**2,947** canonical pairs, **1,317** monomers, 11,768 usable both-direction RR.

### Feasibility statistics (Phase 2)

| metric | value |
|---|---:|
| unique (deduplicated) measurements | **5,742** |
| unique monomer pairs | 2,633 |
| unique monomers | 1,218 |
| condition-complete | **3,225** (56.2 %) |
| quality A / B | 211 / 3,115 |
| multi-source duplicate groups | 302 |
| P2 groups with ≥3 temperatures | 101 (≥4: 41) |
| P2 median temperature span | 45 K |

**P1 Go/No-Go: GO** (≥1500 unique & ≥800 condition-complete).
**P2 feasability: conditional Go** (101 strict multi-T groups).

### P1 benchmark (3,225 condition-complete rows)

Out-of-fold performance of the **best model per split** (5-fold;
`both-monomers-out` test sets are small by construction). Metrics include
calibration (NLL, coverage95) and the Mayo-Lewis `E_curve` (F1-composition
curve error). Full grid (Ridge/RF/XGBoost/MLP × struct/cond/phys) in
`results/p1_benchmark_results.csv`; GNN in `results/p1_gnn_results.csv`.

| target | split | features | model | MAE (ln r) | R² | NLL | cover95 | E_curve |
|---|---|---|---:|---:|---:|---:|---:|---:|
| ln r1 | random | struct+cond | RF | **0.864** | 0.386 | 1.696 | 0.934 | 0.093 |
| ln r1 | pair | struct | RF | 1.061 | 0.220 | 1.816 | 0.940 | 0.120 |
| ln r1 | one-monomer-out | struct+cond+phys | XGB | 1.159 | 0.111 | 1.878 | 0.938 | 0.130 |
| ln r1 | both-monomers-out | struct | RF | 1.202 | 0.076 | 1.912 | 0.931 | 0.136 |
| ln r2 | random | struct+cond | RF | **0.883** | 0.422 | 1.717 | 0.930 | 0.093 |
| ln r2 | pair | struct+cond | RF | 1.059 | 0.301 | 1.812 | 0.934 | 0.116 |
| ln r2 | one-monomer-out | struct+cond+phys | RF | 1.165 | 0.185 | 1.889 | 0.928 | 0.129 |
| ln r2 | both-monomers-out | struct+cond | RF | 1.280 | 0.155 | 1.974 | 0.932 | 0.138 |

**Model comparison (MAE, ln r1):**

| split | Ridge | RF | XGB | MLP | GIN |
|---|---:|---:|---:|---:|---:|
| random | 1.051 | **0.864** | 0.916 | 0.989 | 1.020 |
| pair | 1.252 | **1.061** | 1.076 | 1.164 | 1.079 |
| one-monomer-out | 1.260 | 1.159 | **1.158** | 1.249 | 1.184 |
| both-monomers-out | 1.261 | 1.202 | 1.218 | 1.298 | **1.160** |

**Key findings (P1):**
1. **In-distribution** RR is well predicted from structure (MAE ≈ 0.86–0.88 in
   `ln r`, R² ≈ 0.39–0.42); adding conditions (T, solvent) gives the single
   largest gain in-distribution.
2. **OOD generalisation is the hard problem.** Splitting by pair, one monomer,
   or both monomers roughly doubles the error (MAE 1.06 → 1.28, R² → 0.08–0.30).
   Condition and physics features give little/no OOD gain; error grows
   monotonically with OOD severity (random < pair < one < both).
3. **No model beats the simple fingerprint baseline.** RF/XGB (ECFP4) win almost
   everywhere; the GIN graph network only edges out on the hardest
   `both-monomers-out` split (1.16 vs 1.20) at the cost of in-distribution
   accuracy (1.02 vs 0.86); MLP is uniformly worse.
4. **Calibration tracks OOD degradation.** NLL rises 1.70 → 1.91 and `E_curve`
   0.093 → 0.136 from random to OOD splits, while coverage95 stays ~0.93
   (near-Gaussian, slightly under-confident).
5. **R² is capped by the feature representation, not data noise (corrected).**
   Three independent "clean-the-data" routes all fail to raise R²: (a) pair-level
   aggregation (0.364, fewer rows), (b) noise cleaning — only within-pair MAD
   outlier removal (5 %, genuine unit/transcription errors) gains +0.008 on a
   fixed test while broad quality/uncertainty filtering *hurts* (quality-D
   removal 0.39 → 0.31), and (c) a **primary-source clean subset**
   (`src/p1_clean_subset.py`): `primary` (3494) R² 0.349, `+modern-fit` (479)
   0.313, `+uncertainty` (291) 0.269 — all *below* the full 0.396. The ECFP
   structure model reaches only R² 0.28–0.43, far under the pair-identity
   ceiling (~0.62); the >0.5 gap is **between-pair chemistry** (Q-e polarity,
    sterics, frontier orbitals) that fingerprints cannot express. → Next lever:
    CopDDB Δ-barrier / QM features or Δ-learning residual. Details:
    `src/p1_noise_cleaning.py`, `src/p1_clean_subset.py`,
    `projects/P1_rr_benchmark/p1_noise_pairlevel.json`,
    `projects/P1_rr_benchmark/clean_subset.json`.
 6. **But CopDDB first-principles Δ-barriers do not break the ceiling either**
    (`src/p1_copddb_ceiling.py`): on the 601-measurement / 91-pair / 30-monomer
    CopDDB-covered subset, the Δ-barrier is essentially uncorrelated with
    experimental `ln r` (Pearson +0.18/+0.15, *wrong sign* vs theory) and the
    physical `ln r ≈ −ΔΔE‡/RT` fit gives slope +39 vs the theoretical −1059;
    adding Δ-barrier to ECFP does not raise R² (0.226 → 0.223). → The off-the-shelf
    vacuum small-model DFT barrier proxy carries no experimental signal; breaking
    the ~0.4 ceiling needs genuinely high-quality QM, not existing descriptors.
    Details: `projects/P1_rr_benchmark/copddb_ceiling.json`.

*(See also `projects/P1_rr_benchmark/P1_report.md`.)*

---

## P2 — Temperature-dependent RR thermodynamics

Physics model per (pair, solvent-family, source) series:
`ln r = -ΔΔH‡/(R·T) + ΔΔS‡/R` (linear fit of `ln r` vs `1/T`).

| metric | value |
|---|---:|
| loose series (pair+solvent) | 178 (67 pairs) |
| **strict series (pair+solvent+source)** | **172 (66 pairs)** |
| series with significant T-slope (p<0.05) | 92 (51.7 %) |
| median fit R² | 0.401 |
| median temperature span | 50 K |
| median ΔΔH‡ | 0.0 kJ/mol (IQR −7.4 … +8.5) |
| median ΔΔS‡ | −2.4 J/mol/K |
| enthalpy–entropy compensation slope | **326 K** (r² = 0.99) |

**P2 Go/No-Go: STRONG GO** (172 strict series ≥ 150).

### Key results

1. **Half of all series show a statistically significant temperature response**
   (`ln r` changes with `T`), with a typical |ΔΔH‡| of ~5–8 kJ/mol.
2. **Strong enthalpy–entropy compensation** with `T_comp ≈ 326 K`, i.e. near
   standard polymerization temperatures — **but this is a fit-error artefact.**
   `src/p2_compensation.py` shows the centred regression (which decorrelates the
   OLS slope/intercept) collapses the correlation from r²=0.986 to **r²=0.033**
   (Pearson −0.18), and the apparent compensation temperature (326 K) tracks the
   dataset mean temperature (320 K) — the classic signature of the artefact
   (Krug–Hunter–Grieger). So there is **no physical ΔΔH‡–ΔΔS‡ compensation**
   here; the raw enthalpy variation is real, not cancelled by entropy.
3. **Negative result:** ΔΔH‡ and ΔΔS‡ are **not predictable from monomer
   structure** (leave-pair-out R² ≤ 0 for both Ridge and RF; sign of ΔΔH‡
   predicted at chance, 0.50). Even the 89 significant series give R² = −0.5.
   → Temperature response is effectively *irreducible*; it must be **measured**,
   which directly motivates P3's active experimental design.

Scripts: `src/p2_temperature.py`, `src/p2_thermo_model.py`,
`src/p2_compensation.py`; outputs `projects/P2_temperature_rr/temperature_series_fits.csv`,
`p2_summary.json`, `p2_thermo_model.json`, `p2_compensation_check.json`.
Report: `projects/P2_temperature_rr/P2_report.md`.

## P3 — ML prior + Bayesian / active experimental design

Problem: given a fixed set of **unseen chemistry** (25 % of monomer pairs held
out, 831 test rows), choose which copolymerisation experiments to run under a
budget to minimise prediction error on the held-out pairs. Candidate pool = the
other 75 % of pairs; **shared** cold start = 120 experiments; 12 acquisition
rounds × 30 experiments (final budget 480); model = 5-seed RandomForest
ensemble (predictive mean + std); metric = MAE(ln r1) on held-out pairs.
Reported as **mean ± std over 3 seeds**, plus normalised improvement (share of
the *greedy-oracle* ceiling that a strategy recovers).

| acquisition | MAE @120 | MAE @480 | improvement | norm. improvement |
|---|---:|---:|---:|---:|
| random | 1.255 | 1.241 | 0.014 | 0.18 |
| uncertainty | 1.255 | 1.211 | 0.045 | 0.56 |
| diversity (max-min ECFP) | 1.255 | 1.215 | 0.040 | 0.51 |
| **uncertainty + diversity** | 1.255 | **1.185** | **0.071** | **0.89** |
| physics-prior (CopDDB) | 1.255 | 1.254 | 0.001 | 0.02 |
| greedy-oracle (upper bound) | 1.255 | 1.176 | 0.079 | 1.00 |

**Key result:** the combined **uncertainty + diversity** acquisition is the
winner, reaching **89 % of the greedy-oracle ceiling** (MAE 1.255 → 1.185) —
roughly twice the normalised gain of uncertainty or diversity alone, and ~5×
random. NLL and coverage track the improvement (best strategy NLL 1.90 vs
1.95 for random). The CopDDB physics-prior is ineffective, but this is **not**
evidence that physics is uninformative: CopDDB descriptors cover only **5 % of
candidate pairs (58/1161, 27 monomers)**, so the strategy collapses to picking
from a narrow chemistry subset (a diversity confound, not a physics signal).
This is the quantitative case for an *uncertainty-driven,
diversity-regularised* experimental loop when extending RR coverage to new
monomers — and it directly complements the P2 finding that the temperature
response is irreducible and therefore must be measured.

Script: `src/p3_active_design.py`; outputs
`projects/P3_active_rr_design/active_learning_curve.csv`, `p3_summary.json`.
Report: `projects/P3_active_rr_design/P3_report.md`.

---

## P4 — RR → Sequence → Tg multi-scale Δ-learning

Components (all runnable):

1. **Mayo-Lewis sequence engine** (`src/p4_sequence.py`). Maps `(r1, r2, f1)` to
   instantaneous composition `F1` and terminal-model sequence descriptors
   (P11/P12/P21/P22, dyad fractions, run lengths `L1/L2`, blockiness,
   alternation, sequence entropy). **Validated against Monte-Carlo chains**
   (e.g. r1=1, r2=2, f1=0.5: analytic `F1=0.400` vs MC `0.3997`; `P11` 0.500 vs
   0.502; `P22` 0.667 vs 0.668). Applied to the **2,633 real RR pairs** →
   50,046 feed-composition rows (`pair_sequence_descriptors.parquet`).

2. **Real Tg models from repeat-unit structure** (`src/p4_tg_model.py`).
   PolyInfo repeat-unit SMILES → ECFP4 + RDKit descriptors → Ridge/RF, with a
   chemistry-aware scaffold-cluster split:

   | dataset | n | split | RF MAE (°C) | RF R² |
   |---|---:|---|---:|---:|
   | homopolymers | 8,000 | random | 21.2 | 0.887 |
   | homopolymers | 8,000 | scaffold (OOD) | 34.4 | 0.752 |
   | alt-copolymers | 34 matched | random | 42.7 | 0.414 |

   The homopolymer Tg model generalises (scaffold R² 0.75); the alternating
   copolymer subset matched to two homopolymer components is only **34 rows**,
   because PolyInfo `-alt-` entries are overwhelmingly step-growth
   (polyamide/polyimide/polyurethane) rather than radical copolymers.

3. **Δ-learning pipeline** (`src/p4_delta_learning.py`). Real `r1, r2`
   (RRMaster) + real homopolymer `Tg1, Tg2` (PolyInfo, **254 matched radical
   pairs**, 1,270 composition points) → Mayo-Lewis sequence → Gordon-Taylor/Fox
   baseline; an ML model then learns the sequence correction. Because
   experimental radical-copolymer Tg with composition is not yet in the data,
   the *target* uses a documented dyad-interaction model (γ·(blockiness−1)·ΔT,
   σ=3 °C) — a controlled simulation that validates the methodology:

   | model | MAE (°C) | R² (residual) |
   |---|---:|---:|
   | Fox / Gordon-Taylor baseline | 4.19 | — |
   | + Δ-learning (Ridge) | 3.24 | 0.43 |
   | + Δ-learning (RF) | **3.02** | **0.51** |

   Top residual features are `alternation`/`blockiness` and the `P21/P22`/
   `dyad_22` sequence terms — exactly the sequence dependence built into the
   target. This remains a controlled *sanity check* that the Δ-learning
   machinery can recover a known sequence signal when one exists.

4. **Real radical-copolymer Tg Δ-learning** (`src/p4_real_sequence_tg.py`) —
   the definitive test. Real copolymer Tg data
   (`tgvalue-0416.xlsx` acrylate/methacrylate/St/AN/VC, 351 rows;
   `cleaned-2006-styrenic-co48.xlsx` styrene copolymers, 48 rows) were matched
   to real reactivity ratios from RRMaster → **210 measurements, 40 monomer
   pairs** with `r1·r2 ∈ [0.02, 2.85]`. Leave-one-pair-out CV on the Fox
   residual:

   | model / features | MAE (K) | R² |
   |---|---:|---:|
   | Fox (mole fraction) baseline | 8.60 | 0.925 |
   | Fox (weight) baseline | 9.28 | 0.920 |
   | Δ (composition only) RF | 10.0 | 0.902 |
   | Δ (+ reactivity ratios) RF | 9.90 | 0.896 |
   | Δ (+ sequence descriptors) RF | 10.36 | 0.900 |
   | Δ (+ rr + structure ECFP) RF | **8.14** | **0.923** |

   **Headline result:** the composition (Fox) rule already captures **~92 %** of
   radical-copolymer Tg variance. Reactivity-ratio-derived sequence
   information (`r1·r2` alternation tendency, blockiness, dyad fractions)
   carries **no additional predictive power** — the row-level correlation with
   the Fox residual is `|r| ≤ 0.21`, and the pair-level
   `r1·r2 ↔ mean-residual` correlation is **0.06**. The residual is dominated by
   *monomer-specific chemistry* (e.g. H-bonding in styrene–acrylamide,
   `+25 K`), which is why structure features (ECFP) recover a little more than
   sequence features do. Within the near-random regime covered by common
   radical pairs (`r1·r2` mostly 0.2–0.9), **Tg is essentially
   composition-determined and sequence effects are second-order.**

   This is a useful negative/quantitative result: it bounds how much
   "sequence engineering" can move copolymer Tg for ordinary random
   copolymers, and pins the Fox rule as the correct zero-order baseline.

   **Limitations / next steps:** only 2 of 40 pairs are strongly alternating
   (`r1·r2 < 0.1`) and 2 blocky (`r1·r2 > 1.5`), and those are confounded with
   monomer-specific chemistry (acrylonitrile, acrylamide H-bonding); the
   reported composition columns are copolymer (not feed) composition, so the
   sequence descriptors assume the near-random approximation `F1 ≈ fA`. A
   definitive test of the sequence term needs (a) strongly alternating
   pairs whose comonomer is otherwise inert (e.g. St/fumaronitrile,
   St/maleimide) and (b) feed-vs-copolymer composition explicitly resolved.

5. **End-to-end validation on strong-sequence systems**
   (`src/p4_sequence_tg_validation.py`). The `DeltaLearningCopolymerTg` data
   (item 4) is dominated by near-random pairs, so it could only show a *null*
   sequence effect. The sibling `SequenceTg` repo supplies the missing
   strong-sequence radical systems (Johnston 1973/1976; Barton 1970; NMR-dyad
   radical copolymers). Three linked validations:

   **A. RR → sequence** (Mayo-Lewis vs NMR dyads, 6 systems). The terminal-model
   engine reproduces measured dyad fractions with **MAE 0.055** (hetero-dyad
   MAE 0.043); near-perfect for AN-EA, S-Acrylamide, NVP-CEVE, VDF-HFP, with
   the one systematic miss being MMA-St "random" (RAFT-made, truly random,
   while `r_MMA·r_S=0.24` implies mild alternation).

   **B. RR → sequence → Tg** (Barton dyad-additivity vs Fox, Johnston 5 systems,
   31 samples). The dyad model (dyads from Mayo-Lewis where `r` is known, else
   ideal-random) **halves the error**: Fox MAE 15.3 K → dyad MAE 7.8 K,
   R² 0.754 → 0.881. Per system: S-AN 6.6→0.9 K, AN-MMA 13.6→1.6 K, aMS-AN
   25.0→4.5 K (the only worse case is MMA-VCl, whose alternating dyad is 41 K
   *below* the homopolymer mean, so the single-Tg_AB model over-corrects).

   **C. Barton sign rule.** The direction of the alternating-dyad shift is
   **system-specific**: Tg_AB is *below* the homopolymer mean for 4/5 systems
   (alternation lowers Tg, up to −41 K) but *above* it for S-AN (+9 K). This
   explains both why a universal "alternation → higher Tg" heuristic fails and
   why item 4's near-random pairs showed no sequence signal (their dyads ≈
   ideal-random, so Fox ≈ Barton).

   **Headline (mechanism-level):** the RR → sequence → Tg chain is
   *mechanistically validated* — RR predicts NMR dyads (5 % MAE), and those
   dyads drive a Barton model that halves Fox error on the Johnston systems.
   This is a **within-system descriptive comparison** (the `Tg_AB` parameter is
   reused from the same systems), *not* an independent cross-system prediction.

6. **Supplementary matched-case validation + novelty pressure test** — the
   sharpened, defensible P4 mainline. Full record:
   `projects/P4_sequence_tg/P4_supplementary_experiments.md` (S0–S12),
   `sequence_tg_literature_case_audit.md`, `p4_two_stage_evidence_matrix.csv`.

   * **Novelty boundary.** "Sequence affects Tg" and "RR-derived dyads correct
     Tg" are **not new** (Barton 1970 — 11 systems; SAN 2000; S/EA 1993;
     AN/EA 1992; MMA/EA 2001 report conversion-dependent NMR-vs-RR dyads).
     The defensible claim is the **two-stage error-propagation + applicability
     boundary**: (i) can RR + available feed/composition predict NMR dyads, and
     (ii) do the resulting sequence errors × dyad-Tg leverage yield a
     *transferable* Fox improvement under leave-source / leave-pair splits?
   * **Matched cases (mixed, no blind cross-pair success yet):** SAN
     microemulsion (5 Tg points) RR→heterodyad MAE 0.082, leave-one-out Fox
     3.22 °C vs RR-dyad 1.07 °C (**positive**); S/MA bulk (3 Tg points)
     RR-dyad Tg *worse* than Fox (3.70 vs 2.78 K, **negative boundary**); E/V
     2026 cross-sample (Sample 2→3) dyad MAE 0.159 pp / triad 0.725 pp but fixes
     the target final composition (**conditional**, not blind).
   * **External consistency:** PLGA 2024 (12 samples, 4 fixed compositions)
     shows within-composition alternation↔Tg association (Holm-corrected
     p = 0.042) but is ring-opening, not radical.
   * **Regime map** (`p4_sequence_regime_map.csv`): sequence leverage
     `ΔTg_seq ≈ ΔH × ΔT_AB` (departure from Bernoulli × dyad-Tg interaction)
     bounds the effect to −9.5…+4.5 K over the Johnston 5 systems — only systems
     with both large `ΔH` *and* large `ΔT_AB` can expect a visible sequence
     correction.
   * **Positioning:** a bounded matched-case synthesis / method validation, not
     a generalisation claim (no cross-pair, no-target-leakage end-to-end
     success yet; minimum bar in S12).

Scripts: `src/p4_sequence.py`, `src/p4_tg_data.py`, `src/p4_tg_model.py`,
`src/p4_delta_learning.py`, `src/p4_real_sequence_tg.py`,
`src/p4_sequence_tg_validation.py`, `src/p4_supplementary_plga.py`,
`src/p4_supplementary_feed_sensitivity.py`, `src/p4_supplementary_rr_spread.py`,
`src/p4_condition_and_source_cv.py`, `src/p4_literature_SAN_case.py`,
`src/p4_literature_SMA_case.py`, `src/p4_ev_c2_sequence_holdout.py`,
`src/p4_sequence_regime_map.py`; outputs under `projects/P4_sequence_tg/`.
Report: `projects/P4_sequence_tg/P4_report.md`.

---

## Shared data-reuse assets

| asset | path |
|---|---|
| Monomer registry | `metadata/monomer_registry.csv` |
| Solvent registry | `metadata/solvent_registry.csv` |
| Source lineage | `metadata/source_lineage.csv` |
| Duplicate groups | `data/provenance/duplicate_groups.csv` |
| RRMaster v1 | `data/standardized/rrmaster_v1.parquet` |

## Change log

* **2026-10-06** — P4 supplementary matched-case validation + novelty pressure
  test (S0–S12): literature audit shows "sequence affects Tg" and "RR-derived
  dyads correct Tg" are **not new** (Barton 1970, SAN 2000, S/EA 1993, AN/EA
  1992, MMA/EA 2001); reframed the P4 mainline to a **two-stage
  error-propagation + applicability boundary** study. Added matched cases SAN
  (positive), S/MA (negative), E/V 2026 (conditional cross-sample), PLGA 2024
  (external), plus a sequence-leverage regime map (`ΔTg_seq ≈ ΔH × ΔT_AB`).
  Updated `P4_report.md` and README P4 §6.
* **2026-10-05** — scaffolded project; Stage A–G pipeline → RRMaster v1;
  Phase-2 feasibility (P1 GO, P2 conditional GO); P1 benchmark smoke run.
* **2026-10-05** — completed P1 full benchmark (random/pair/one- & both-monomer
  out, RF/XGB/Ridge, condition + CopDDB physics ablation). Key result: OOD
  error ~2× in-distribution and dominated by chemistry coverage.
* **2026-10-05** — completed P2 (172 strict multi-T series, **strong Go**;
  structure→ΔH‡/ΔS‡ is a negative result). Completed P3 (uncertainty-driven
  active design wins). Completed P4 framework+demo (MC-validated sequence
  engine; homopolymer Tg R² 0.75 scaffold; Δ-learning 4.19→3.02 °C).
* **2026-10-05** — P4 real-data upgrade: matched 210 real radical-copolymer Tg
  measurements to real reactivity ratios (40 pairs). **Key result: Fox rule
  already explains ~92 % of Tg variance; reactivity-ratio sequence information
  adds nothing (r1·r2 ↔ residual corr ≈ 0.06) — Tg is composition-determined
  for ordinary random copolymers.**
* **2026-10-05** — P4 definitive validation on strong-sequence systems
  (`SequenceTg` Johnston/Barton + NMR dyads): RR → sequence reproduces NMR dyads
  (MAE 0.055); Barton dyad model **halves Fox error** (15.3→7.8 K, R² 0.75→0.88);
  the alternating-dyad shift is **system-specific** (Tg_AB below homopolymer
  mean in 4/5 systems, above in S-AN).
* **2026-10-05** — P1 finalised: GNN (GIN) + MLP baselines, calibration (NLL /
  coverage) + Mayo-Lewis E_curve metrics, OOF preds saved; wrote `P1_report.md`.
  **No model beats ECFP+RF/XGB; OOD error dominated by chemistry coverage.**
* **2026-10-05** — P2 compensation verified as a **fit-error artefact**
  (centred r² 0.986 → 0.033); wrote `P2_report.md`. P3 rebuilt as a strict
  3-seed offline benchmark (shared cold start, calibration, greedy-oracle
  ceiling): **uncertainty+diversity wins (89 % of oracle)**, CopDDB physics
  prior useless; wrote `P3_report.md`.
* **2026-10-05** — P1 noise-floor correction + clean-subset: three "clean the
  data" routes (replicate count, MAD outliers, primary-source subset) all fail
  to raise R² → binding constraint is the **feature representation** (ECFP gap
  for between-pair chemistry), not measurement noise. Wrote `p1_clean_subset.py`;
  updated `P1_report.md` §3.7 + conclusions #5/#6 and README key-finding #5.
* **2026-10-05** — P1 ceiling verification (`src/p1_copddb_ceiling.py`): CopDDB
  first-principles Δ-barriers are uncorrelated with experimental `ln r` (Pearson
  +0.18, wrong sign; slope +39 vs theory −1059) and add nothing to ECFP → the
  ~0.4 R² ceiling is a representation limit not solved by off-the-shelf DFT
  barriers. Updated `P1_report.md` §3.8 + conclusion #7, README key-finding #6.
* **2026-10-05** — Coxwell curated subset (`src/p1_clean_subset.py` now reads the
  raw xlsx curation flags): removing the Coxwell authors' own secondary rows
  (polymer handbook / Young compilation / "need to delete") yields 220 rows but
  R² **falls** (0.366→0.316 ln_r1, 0.445→0.396 ln_r2) — a fourth confirmation
  that data cleaning is not the bottleneck.
* **2026-10-05** — P1 noise-floor + cleaning analysis
  (`src/p1_noise_cleaning.py`): irreducible within-pair scatter caps R² at
  ~0.62; pair-level aggregation does not raise R²; only within-pair MAD-outlier
  removal helps (marginally), broad quality/uncertainty filtering hurts.
