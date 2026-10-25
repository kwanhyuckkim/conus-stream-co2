"""
Parameter sweep for a single basin (Matt's proposed experiment).

Goal
  Before committing to LHS/Optuna/JAX, test two things on ONE basin:
    (1) linearity: do the basin medians move (nearly) linearly with each knob?
    (2) influence: how far does each knob actually move the medians across its
        plausible range? (Matt suspected Fwc may barely move the basin scale.)

IMPORTANT wording
  We report INFLUENCE (marginal one-at-a-time sensitivity), NOT identifiability.
  A knob can be influential yet non-identifiable if its effect is collinear with
  another knob's; that only shows up in the 2-D (joint) sweep, not in a 1-D slice.

This module does NOT change the model. It only calls calibrate.evaluate() many times.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import calibrate
from .constants import (FWC_LAKE_LOWER_PHOTO, FWC_LAKE_UPPER_PHOTO, GREAT_LAKES_HUC4,
                        LOWER_BOUNDS, PARAM_NAMES, PHOTOSYNTHESIS_HUC2, UPPER_BOUNDS)

# Below this span [ppm] a knob is treated as having no effect (flat), so linearity is undefined.
INFLUENCE_FLOOR_PPM = 1.0

# Water-body a knob is naturally judged on (its own domain).
OWN_RESPONSE = {"Cbz_riv": "river_median", "Fwc_riv": "river_median",
                "Cbz_lake": "lake_median", "Fwc_lake": "lake_median"}
CROSS_RESPONSE = {"river_median": "lake_median", "lake_median": "river_median"}


def _norm_huc4(huc4):
    """Coerce HUC4 to str (NHD codes have leading zeros) and reject Great Lakes basins."""
    if huc4 is None:
        return None
    if isinstance(huc4, (int, float, np.integer, np.floating)):
        huc4 = str(int(huc4))          # numeric HUC4 (418, 418.0) -> '418'
    huc4 = str(huc4)
    if huc4.isdigit():
        huc4 = huc4.zfill(4)   # recover a leading zero an int HUC4 dropped (418 -> '0418')
    if huc4 in GREAT_LAKES_HUC4:
        raise ValueError(
            f"HUC4 {huc4!r} is a Great Lakes basin. The R pipeline excludes these (returns all-NA) "
            f"before running the model (calibrate_model.R:42), so a sweep is not defined for them."
        )
    return huc4


def knob_bounds(knob_index: int, huc4: str | None = None) -> tuple[float, float]:
    """Return (lower, upper) for one knob, applying the HUC2 04/16 photosynthesis override."""
    huc4 = _norm_huc4(huc4)
    lo, hi = LOWER_BOUNDS[knob_index], UPPER_BOUNDS[knob_index]
    if knob_index == 3 and huc4 is not None and huc4[:2] in PHOTOSYNTHESIS_HUC2:
        lo, hi = FWC_LAKE_LOWER_PHOTO, FWC_LAKE_UPPER_PHOTO   # net lake photosynthesis
    return lo, hi


def default_baseline(huc4: str | None = None) -> tuple[float, float, float, float]:
    """Neutral baseline = midpoint of each knob's range. NOTE: this is an UNCALIBRATED operating
    point, so absolute pCO2 along a slice can be large/nonphysical; read the SHAPE, not the value."""
    huc4 = _norm_huc4(huc4)
    return tuple(0.5 * (lo + hi) for lo, hi in (knob_bounds(k, huc4) for k in range(4)))


def sweep_1d(basin, riverCO2, lakeCO2, baseline, knob_index, n=25, huc4=None, upstream_df=None):
    """Vary ONE knob across its bounds (others held at baseline). Return a tidy DataFrame.

    Columns: knob, value, river_median, lake_median, cost, fitness
    upstream_df: cross-basin CO2 imports (usually None for a single/headwater basin)."""
    huc4 = _norm_huc4(huc4)
    lo, hi = knob_bounds(knob_index, huc4)
    rows = []
    for v in np.linspace(lo, hi, n):
        p = list(baseline)
        p[knob_index] = float(v)
        ev = calibrate.evaluate(basin, tuple(p), riverCO2=riverCO2, lakeCO2=lakeCO2,
                                huc4=huc4, upstream_df=upstream_df)
        rows.append({"knob": PARAM_NAMES[knob_index], "value": float(v),
                     "river_median": ev.river_median, "lake_median": ev.lake_median,
                     "cost": ev.cost, "fitness": ev.fitness})
    return pd.DataFrame(rows)


def sweep_all_1d(basin, riverCO2, lakeCO2, baseline=None, n=25, huc4=None, upstream_df=None):
    """Run a 1-D (one-at-a-time) sweep for all four knobs. Returns {knob_name: DataFrame}."""
    huc4 = _norm_huc4(huc4)
    if baseline is None:
        baseline = default_baseline(huc4)
    return {PARAM_NAMES[k]: sweep_1d(basin, riverCO2, lakeCO2, baseline, k, n=n,
                                     huc4=huc4, upstream_df=upstream_df) for k in range(4)}


def influence_ppm(df: pd.DataFrame, response_col: str) -> float:
    """Marginal sensitivity: how far `response_col` moves across the knob's full range [ppm].

    IMPORTANT: this ppm span is a VALUE conditional on where the other three knobs are held (the
    midpoint baseline in default_baseline), not a baseline-invariant quantity. With interactions a
    knob can look more/less influential purely because of the hold-point. The 2-D sweep is the
    authority on coupling and cross-knob ranking; treat this number as "at this baseline"."""
    y = df[response_col].to_numpy(dtype=float)
    y = y[np.isfinite(y)]
    # need >=2 finite points for a span; 1 finite point (rest NaN) is undefined, not "flat".
    return float(y.max() - y.min()) if len(y) >= 2 else np.nan


def linearity_r2(df, response_col, influence_floor=INFLUENCE_FLOOR_PPM):
    """R^2 of a straight-line fit of `response_col` vs the knob value. 1.0 = perfectly linear.

    Returns NaN when the knob barely moves the response (span < influence_floor): a dead/flat knob
    is NOT 'linear', it is uninformative, and scoring it 1.0 would inflate an 'everything is linear'
    read-out (exactly the Fwc case Matt flagged).

    NOTE: response is a MEDIAN (an order statistic), so this tests the basin-median summary, not the
    per-reach physics."""
    x = df["value"].to_numpy(dtype=float)
    y = df[response_col].to_numpy(dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return np.nan
    span = y.max() - y.min()
    if span < influence_floor:            # flat / dead knob -> linearity undefined
        return np.nan
    slope, intercept = np.polyfit(x, y, 1)
    y_hat = slope * x + intercept
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def summarize_1d(sweeps: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """One row per knob: own-domain influence, cross-domain influence (exposes river/lake coupling),
    and linearity R^2 (NaN if the knob is flat)."""
    rows = []
    for name, df in sweeps.items():
        own = OWN_RESPONSE[name]
        cross = CROSS_RESPONSE[own]
        infl = influence_ppm(df, own)
        r2 = linearity_r2(df, own)                       # compute once
        rows.append({
            "knob": name, "judged_on": own,
            "influence_ppm": round(infl, 3) if np.isfinite(infl) else np.nan,
            "cross_influence_ppm": round(influence_ppm(df, cross), 3),  # coupling
            "linearity_R2": (np.nan if not np.isfinite(r2) else round(r2, 5)),
            "flat": bool(np.isfinite(infl) and infl < INFLUENCE_FLOOR_PPM),
        })
    return pd.DataFrame(rows)


def sweep_2d(basin, riverCO2, lakeCO2, baseline, knob_i, knob_j, n=20, huc4=None, upstream_df=None):
    """Grid two DISTINCT knobs (others at baseline); return long-form cost/median surface."""
    if knob_i == knob_j:
        raise ValueError("sweep_2d requires two distinct knobs (knob_i != knob_j).")
    huc4 = _norm_huc4(huc4)
    gi = np.linspace(*knob_bounds(knob_i, huc4), n)
    gj = np.linspace(*knob_bounds(knob_j, huc4), n)
    rows = []
    for vi in gi:
        for vj in gj:
            p = list(baseline)
            p[knob_i] = float(vi)
            p[knob_j] = float(vj)
            ev = calibrate.evaluate(basin, tuple(p), riverCO2=riverCO2, lakeCO2=lakeCO2,
                                    huc4=huc4, upstream_df=upstream_df)
            rows.append({PARAM_NAMES[knob_i]: float(vi), PARAM_NAMES[knob_j]: float(vj),
                         "cost": ev.cost, "fitness": ev.fitness,
                         "river_median": ev.river_median, "lake_median": ev.lake_median,
                         # per-target error DIAGNOSTICS only. NOTE: the off-body term is NOT constant
                         # on a same-body slice -- lake knobs propagate into downstream river reaches
                         # via routing (and vice versa), so judge the behavioral set on the full 'cost'
                         # (= river_abs_err + lake_abs_err), the same 2-target objective the GA uses.
                         "river_abs_err": abs(ev.river_median - riverCO2),
                         "lake_abs_err": abs(ev.lake_median - lakeCO2)})
    return pd.DataFrame(rows)


def behavioral_set(surface_2d: pd.DataFrame, abs_tol_ppm: float = 50.0, cost_col: str = "cost"):
    """Keep points within an ABSOLUTE ppm budget of the slice's best cost (GLUE-style).

    Absolute (best + abs_tol_ppm), NOT relative (best*(1+tol)): the relative form degenerates -- it
    fills the whole slice when targets are unreachable (large best) and collapses to the exact optimum
    when best is 0. Note: `best` is this slice's minimum, which need not reach the global optimum, so
    print it."""
    best = float(surface_2d[cost_col].min())
    if not np.isfinite(best):
        raise ValueError(f"behavioral_set: cost column {cost_col!r} is empty or all-NaN "
                         "(both water-bodies absent, or medians NaN-poisoned).")
    # The slice optimum (cost==best) is ALWAYS retained, so this set is never empty. A size-1 result
    # means abs_tol_ppm is below the cost change per grid step, NOT that targets are unreachable.
    keep = surface_2d[surface_2d[cost_col] <= best + abs_tol_ppm].copy()
    keep.attrs["slice_best_cost"] = best
    return keep
