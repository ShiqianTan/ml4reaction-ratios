"""Conditional Sample 2 -> Sample 3 sequence prediction for the 2026 E/V case.

This test uses Sample 2's C2 reactivity ratios and Sample 3's reported polymer
composition. It predicts sequence conditional on composition; it is not a blind
feed-to-polymer prediction because the supporting-information feed conditions
are not yet available.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "projects" / "P4_sequence_tg" / "ev_2026_matched_case.csv"
OUT_CSV = ROOT / "projects" / "P4_sequence_tg" / "ev_2026_c2_sample2_to_3_sequence_holdout.csv"
OUT_JSON = ROOT / "projects" / "P4_sequence_tg" / "p4_ev_c2_sequence_holdout.json"


def terminal_transition_probabilities(r_e: float, r_v: float, feed_v: float) -> tuple[float, float]:
    """Return q(E->V), q(V->E) for the terminal Mayo-Lewis model."""
    feed_e = 1.0 - feed_v
    q_ev = feed_v / (r_e * feed_e + feed_v)
    q_ve = feed_e / (feed_e + r_v * feed_v)
    return q_ev, q_ve


def feed_v_for_polymer_composition(r_e: float, r_v: float, target_v: float) -> float:
    """Invert the stationary Markov composition relation by bisection."""
    lo, hi = 0.0, 1.0
    for _ in range(100):
        feed_v = (lo + hi) / 2.0
        q_ev, q_ve = terminal_transition_probabilities(r_e, r_v, feed_v)
        polymer_v = q_ev / (q_ev + q_ve)
        if polymer_v < target_v:
            lo = feed_v
        else:
            hi = feed_v
    return (lo + hi) / 2.0


def predict(r_e: float, r_v: float, polymer_v: float) -> tuple[dict[str, float], dict[str, float], float]:
    """Predict stationary terminal-model dyad and triad fractions."""
    feed_v = feed_v_for_polymer_composition(r_e, r_v, polymer_v)
    q_ev, q_ve = terminal_transition_probabilities(r_e, r_v, feed_v)
    q_ee, q_vv = 1.0 - q_ev, 1.0 - q_ve
    polymer_e = 1.0 - polymer_v

    dyads = {
        "EE": polymer_e * q_ee,
        "EV": polymer_e * q_ev + polymer_v * q_ve,
        "VV": polymer_v * q_vv,
    }
    triads = {
        "EEE": polymer_e * q_ee * q_ee,
        "EEV": 2.0 * polymer_e * q_ee * q_ev,
        "VEV": polymer_v * q_ve * q_ev,
        "EVE": polymer_e * q_ev * q_ve,
        "EVV": 2.0 * polymer_v * q_ve * q_vv,
        "VVV": polymer_v * q_vv * q_vv,
    }
    return dyads, triads, feed_v


def mae_rmse(errors: list[float]) -> tuple[float, float]:
    return sum(abs(e) for e in errors) / len(errors), math.sqrt(sum(e * e for e in errors) / len(errors))


def main() -> None:
    with DATA.open(newline="", encoding="utf-8") as stream:
        samples = {row["sample"]: row for row in csv.DictReader(stream)}
    train = samples["2"]
    target = samples["3"]
    r_e, r_v = float(train["rE_reported"]), float(train["rV_reported"])
    polymer_v = float(target["xV_author_text_mol_pct"]) / 100.0
    dyads, triads, effective_feed_v = predict(r_e, r_v, polymer_v)

    observed = {
        "EE": 0.685,
        "EV": 0.265,
        "VV": 0.050,
        "EEE": 0.581,
        "EEV": 0.210,
        "VEV": 0.029,
        "EVE": 0.103,
        "EVV": 0.061,
        "VVV": 0.016,
    }
    rows = []
    for group, predictions in (("dyad", dyads), ("triad", triads)):
        for sequence, prediction in predictions.items():
            measured = observed[sequence]
            rows.append({
                "prediction_group": group,
                "sequence": sequence,
                "observed_sample3_pct": measured * 100.0,
                "predicted_from_sample2_rr_pct": prediction * 100.0,
                "error_percentage_points": (prediction - measured) * 100.0,
            })
    with OUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    dyad_errors = [(dyads[k] - observed[k]) * 100.0 for k in dyads]
    triad_errors = [(triads[k] - observed[k]) * 100.0 for k in triads]
    dyad_mae, dyad_rmse = mae_rmse(dyad_errors)
    triad_mae, triad_rmse = mae_rmse(triad_errors)
    result = {
        "test": "conditional terminal-model prediction of Sample 3 sequence from Sample 2 C2 RR",
        "training_sample": 2,
        "target_sample": 3,
        "input_rr": {"rE": r_e, "rV": r_v},
        "target_polymer_v_mol_pct_reported": polymer_v * 100.0,
        "effective_feed_v_inverted_from_target_composition_not_measured_feed": effective_feed_v * 100.0,
        "predicted_dyads_pct": {k: v * 100.0 for k, v in dyads.items()},
        "observed_dyads_sample3_pct": {k: observed[k] * 100.0 for k in dyads},
        "dyad_mae_percentage_points": dyad_mae,
        "dyad_rmse_percentage_points": dyad_rmse,
        "predicted_triads_pct": {k: v * 100.0 for k, v in triads.items()},
        "observed_triads_sample3_pct": {k: observed[k] * 100.0 for k in triads},
        "triad_mae_percentage_points": triad_mae,
        "triad_rmse_percentage_points": triad_rmse,
        "scope_note": "Sequence-only holdout conditional on reported target polymer composition; not a fully blind feed-to-polymer prediction. RR are rounded table values and are themselves estimated from Sample 2 NMR sequence data.",
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
