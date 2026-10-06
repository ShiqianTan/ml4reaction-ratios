# Sequence–Tg literature case audit

## Scope

This is a targeted audit of papers that jointly report reactivity ratios (RR),
NMR-derived sequence statistics, and glass-transition temperature (Tg). “Matched”
means the variables are reported for the same binary monomer system and relevant
sample series; it does not imply that every individual RR estimate was obtained
independently from the Tg samples. Raw row availability is stated separately.

## Cases located

| Pair / paper | RR | Sequence evidence | Tg evidence | Data status for P4 |
|---|---|---|---|---|
| Ethene–vinylcyclohexane (E/V), Cipullo et al. 2026 | Table 4 reports sequence-derived ratios: Sample 1/C1 `rE=19, rV=0.003`; Sample 2/C2 `14, 0.14`; Sample 3/C2 `16, 0.11` | Table 2 gives NMR dyads/triads for 3 samples and Markov-calculated distributions. Sample 1/3 differ in [EV] (35.1 vs 26.5%) and [VV] (0.5 vs 5.0%); paper text reports xV=17.9% (Sample 1) and 18.1% (Sample 3) | Samples 1 and 3: Tg ≈−23 and −14 °C; Mn ≈268 and 401 kDa; DSC also reports Tm=23.1/44.8 °C and ΔHm=5.7/12.7 J g−1 | **Main article and SI transcribed; conditional Sample 2→3 sequence holdout completed.** Table 2 rounded dyads imply xV=18.05/18.25%, while reported composition values differ slightly; preserve both. RR are derived from measured sequence; SI does not report total liquid/solvent volume or ethylene uptake history, preventing end-to-end feed prediction. [Open-access article](https://doi.org/10.1016/j.polymer.2026.130510) |
| Styrene–acrylonitrile (SAN), microemulsion, Sanghvi et al. 2000 | Present-work system-level values reported as `rS=1.49`, `rA=0.022`; table also compares reported bulk, solution, emulsion and microemulsion values | Table 3 reports 13C-NMR triads; Table 4 reports NMR-derived dyads and sequence lengths | Table 4 has 5 experimental DSC Tg values plus a sixth missing value, and reported Barton predictions | **Six rows transcribed.** RR-to-heterodyad and Fox/NMR/RR sequence analyses reproduced from tables; one pair, one paper, n=5 Tg observations. [Article and tables](https://electronicsandbooks.com/edt/manual/Magazine/E/European%20Polymer%20Journal/2000%20Volume%2036/2000%20Volume%2036/10/2275-2283.pdf); [publisher record](https://doi.org/10.1016/S0014-3057(00)00002-1) |
| Styrene–ethyl acrylate, Brar & Sunita 1993 | EVM-derived ratios reported | 13C-NMR triads/diads and calculated sequence parameters; observed triads compared with RR-based values | Barton Tg reported and compared with DSC | **Confirmed by primary abstract; numerical rows not yet transcribed.** [Wiley record](https://doi.org/10.1002/macp.1993.021940618) |
| Acrylonitrile–methyl acrylate, Brar & Sunita 1992 | EVM-derived terminal/penultimate ratios | NMR-derived triad/diad and sequence statistics | DSC Tg compared with NMR-based Tg | **Confirmed by abstract; numerical rows not yet transcribed.** [Wiley record](https://doi.org/10.1002/pola.1992.080301209) |
| Acrylonitrile–ethyl acrylate, Brar & Sunita 1992 | Nonlinear EVM ratios; terminal and penultimate ratios also calculated from NMR triads | 13C-NMR triads converted to dyads and sequence parameters; observed triads reportedly agree with RR-based estimates | Multiple DSC Tg values reportedly agree with dyad-derived values | **Best next older-paper candidate for a third matched case; abstract confirms all data types, but sample-level numerical tables are not accessible in the sources found so far.** [Publisher record](https://doi.org/10.1016/0014-3057(92)90087-I) |
| Styrene–methyl acrylate, Brar & Sunita 1992 | Kelen–Tudos ratios, `rS=1.20`, `rM=0.18` | 13C-NMR triads, dyads, sequence lengths, and block character | Table 4 reports Barton Tg for 6 samples and DSC Tg for 3 | **Six rows transcribed and calculations reproduced.** Terminal RR at reported feed misses NMR heterodyads (MAE 0.266); Tg results depend on whether composition is measured or predicted. [Primary article PDF and tables](https://nopr.niscpr.res.in/bitstream/123456789/45983/1/IJCA%2031A%2812%29%20903-908.pdf) |
| Methyl methacrylate–ethyl acrylate, Nguyen et al. 2001 | Classical reactivity-ratio dyad estimates compared with 13C-NMR dyads; NMR becomes more accurate in high-conversion samples | Compositional dyads plus dyad-level tacticity | Barton/Johnston Tg calculations compared with DSC | **Primary abstract reports that NMR dyads give more precise Tg calculations than RR-derived dyads, especially at high conversion; adding tacticity slightly improves Tg.** Numerical rows not yet transcribed. [Primary article DOI](https://doi.org/10.1002/pi.694) |

The SAN article makes clear that this broad scientific idea is **not itself new**:
sequence-dependent Tg had already been evaluated using measured NMR dyads and a
Barton relation. A publishable gap is narrower: test when **RR-derived** sequence
statistics (rather than NMR sequence measured on the target sample) recover that
signal and improve a Fox baseline under genuinely held-out chemistry or source.
The 2026 E/V study strengthens this novelty boundary: it reports measurable Tg
differences between comparable-composition samples with different NMR-derived
sequences. Its RR are inferred from the same sequence statistics, so it supports
sequence–Tg physics but does not by itself test external RR-driven Tg prediction.
The E/V article's sample 1 and 3 values are approximately −23 and −14 °C for Tg
and 268 and 401 kDa for Mn, respectively. The paper reports xV=17.9% for Sample 1
and 18.1% for Sample 3. Calculating xV from the rounded Table 2 dyads instead
gives 18.05% and 18.25%; these should not be mistaken for more precise measured
compositions. The paper describes molecular weights as comparable, although the
reported values are not identical, so the pair is suggestive evidence for a
sequence contribution rather than a fully isolated causal estimate. Table 4
reports RR derived from sequence analysis; for the shared C2 catalyst, estimates
from Samples 2 and 3 are similar (`rE=14, rV=0.14` vs `rE=16, rV=0.11`).
An independent RR-to-sequence prediction may be testable by fitting on Sample 2
and predicting Sample 3, but requires the supporting-information feed/condition
details and must not use Sample 3 NMR to refit RR.

## Reanalysis of the SAN microemulsion table

Input values in `san_microemulsion_matched_case.csv` were transcribed from Tables
1 and 4 of Sanghvi et al. The study reports anionic oil-in-water microemulsion
polymerization at 70 °C, 5% SAN monomer, 14% SDS, 81% water, and 0.37 mM KPS.
The polymerization-locus styrene fraction was reconstructed as
`M0S/(M0S+M0A)` from the Table 1 `M0S/M0A` ratio. RR-derived dyads use the
paper’s system-level `rS=1.49`, `rA=0.022` in the terminal Mayo–Lewis model.

| Comparison | n | MAE (°C) | RMSE (°C) | Meaning |
|---|---:|---:|---:|---|
| Fox, using measured copolymer composition and paper homopolymer Tg endpoints | 5 | 3.22 | 3.77 | Composition-only baseline |
| Published Barton predictions | 5 | 0.34 | 0.49 | In-sample values computed from measured NMR dyads; not an independent RR prediction |
| Refit dyad interaction in leave-one-sample-out, measured NMR dyads | 5 | 0.78 | 0.90 | Within-paper Tg holdout, sequence measured directly |
| Refit dyad interaction in leave-one-sample-out, RR-derived dyads | 5 | 1.07 | 1.26 | Within-paper Tg holdout, sequence estimated from RR and locus feed |

At the same time, RR-derived heterodyad fractions had mean absolute deviation
**0.082** from the NMR heterodyad fractions across all six samples (0.077 among
the five Tg-evaluable rows). This shows a feasible example where RR recovers the
measured sequence trend and RR-derived dyads retain a Tg signal beyond Fox in a
small within-system leave-one-out exercise.

### Interpretation boundary

- This is one monomer pair, one source, one polymerization medium, and only five
  Tg-evaluable samples. It is a case study, not evidence of generalization.
- The article’s RR are estimated within the same sample series; the held-out Tg
  is held out from the Tg interaction fit, but not from RR estimation.
- The leave-one-out sequence relation estimates an interaction coefficient from
  the remaining four Tg rows. It therefore does not establish transfer to a new
  pair, chemistry, or study.
- The paper itself already demonstrated NMR-dyad/Barton agreement. P4 novelty
  cannot be “sequence affects Tg”; it must be a rigorous scope/transfer test of
  RR-derived sequence descriptors against Fox, including failures and uncertainty.

Reproduction:

```bash
PYTHONPATH=src python3 src/p4_literature_SAN_case.py
```

Outputs: `p4_san_matched_case.json`,
`san_microemulsion_matched_case_results.csv`, and
`san_microemulsion_tg_loo_predictions.csv` in this directory.

## Second matched case: styrene–methyl acrylate

Brar and Sunita report six NMR dyad rows, Kelen–Tudos terminal ratios
(`rS=1.20`, `rM=0.18`), and Tg values; only the first three have experimental
DSC Tg. Samples were made by bulk polymerization at 70 °C and 5–10% conversion.
The article reports homopolymer Tg values of 100 °C (PS) and 6 °C (PMA), and
compares NMR dyads with reactivity-ratio sequence predictions. Values were
transcribed from Tables 3–4 of the [primary article](https://nopr.niscpr.res.in/bitstream/123456789/45983/1/IJCA%2031A%2812%29%20903-908.pdf).

| Analysis | n | Result | Interpretation |
|---|---:|---:|---|
| RR heterodyad vs NMR at reported feed | 6 | MAE 0.266 | Strong terminal-RR sequence mismatch in this case |
| RR heterodyad vs NMR after inverting feed to match NMR-derived copolymer composition | 6 | MAE 0.233 | Still a large mismatch; composition matching alone does not resolve it |
| Fox using measured polymer composition | 3 | MAE 2.78 K | Baseline |
| NMR dyad Tg interaction, sample-level LOO | 3 | MAE 1.37 K | Small within-paper sample holdout; fit on two other Tg points |
| RR dyad at reported feed, Tg interaction LOO | 3 | MAE 3.70 K | Worse than Fox; composition and sequence errors are both present |
| RR dyad at NMR-matched composition, Tg interaction LOO | 3 | MAE 1.36 K | Similar to NMR-dyad fit, but composition was constrained using target NMR data |

This case adds a useful **failure boundary** to the SAN positive case. Feed-based
terminal RR predicts too much heterodyad and worse Tg error than Fox. Forcing the
RR calculation to reproduce the NMR-inferred polymer composition reduces the
Tg error, but this is a diagnostic sensitivity analysis: it uses target-sample
composition and fits one heterodyad Tg interaction using only two training Tg
observations per fold. At `n=3`, the near-tie between the composition-constrained
RR and measured-NMR sequence fits is not robust evidence that RR can recover the
sequence effect. The composition-constrained RR heterodyads remain far from NMR
(MAE 0.233); the similar Tg error shows this three-point Tg analysis cannot
discriminate sequence-reconstruction quality.

The paired cases therefore support a sharper hypothesis: **RR-derived sequence
information can improve on Fox only when the RR/terminal model and the
polymerization-locus composition jointly reproduce the sample sequence, and the
Tg sequence correction transfers beyond the samples used to fit it.** SAN meets
these criteria descriptively within one study; S/MA does not when evaluated from
the reported feed. Neither case establishes transfer to a held-out monomer pair.

Reproduction:

```bash
PYTHONPATH=src python3 src/p4_literature_SAN_case.py
PYTHONPATH=src python3 src/p4_literature_SMA_case.py
```

S/MA input and outputs: `styrene_methyl_acrylate_matched_case.csv`,
`styrene_methyl_acrylate_matched_case_results.csv`,
`styrene_methyl_acrylate_tg_loo_predictions.csv`, and
`p4_sma_matched_case.json`.

## E/V supporting information and conditional sequence holdout

The local SI DOCX reports Table S1 conditions: all runs used 80 °C and 80 psi
ethylene, with MAO/BHT cocatalyst at 1:0.4. Sample 1 used catalyst C1 (6 nmol
Hf), 1.00 mL VCH, 520 s, and yielded 149 mg; Sample 2 used C2 (6 nmol Hf),
2.00 mL VCH, 1800 s, and yielded 100 mg; Sample 3 used C2 (1 nmol Hf), 0.75 mL
VCH, 2300 s, and yielded 82 mg. The machine-readable transcription is in
`ev_2026_s1_conditions.csv`.

The SI also gives the NMR integration relationships for dyads/triads and the
equations used to derive reactivity ratios from sequence statistics. Thus the
reported RR are sequence-fitted quantities, not an independent measurement
that avoids NMR.

Using Sample 2's C2 ratios (`rE=14, rV=0.14`) and Sample 3's reported polymer
composition (`xV=18.1 mol%`) as a fixed input, the terminal Mayo–Lewis model
predicts Sample 3 dyad fractions with MAE **0.159 percentage points** (RMSE
0.176) and triad fractions with MAE **0.725 points** (RMSE 0.790). Predicted
EV/VV are 26.32/4.94%, versus measured 26.5/5.0%. Sample 3 NMR was not used to
refit RR, making this a useful cross-sample sequence holdout conditional on
composition.

This is not a feed-to-sequence blind prediction. Table S1 omits solvent/total
liquid volume and ethylene uptake (including its time course). Yield and final
composition can retrospectively constrain incorporated monomer mass, but they
do not provide the pre-reaction inputs or the evolving dissolved ethylene
concentration needed to reconstruct the semi-batch feed history. See
`ev_2026_c2_sample2_to_3_sequence_holdout.csv`,
`p4_ev_c2_sequence_holdout.json`, and
`src/p4_ev_c2_sequence_holdout.py` for data and calculation.
