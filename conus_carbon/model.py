"""
CO2_MODEL: single-reach steady-state CO2 transport (the physics core).

Faithful port of `src/model.R` (CO2_MODEL, equation 6).

FLOW IS UNCHANGED
  The equation and every branch match the R code line-by-line. The ONLY change
  is that the caller passes the already-gathered upstream CO2/Q arrays instead of
  making this function scan `toNode_vec == fromNode` itself. That scan was O(N)
  per reach (O(N^2) total); precomputing an upstream-index map upstream makes it
  O(1) lookup here, with identical numerical results.
"""
from __future__ import annotations

import numpy as np


def _weighted_mean_narm(x: np.ndarray, w: np.ndarray) -> float:
    """R's weighted.mean(x, w, na.rm=TRUE): drop entries where x is NaN, then sum(x*w)/sum(w)."""
    mask = ~np.isnan(x)                       # keep non-NaN CO2
    xw = x[mask]
    ww = w[mask]
    denom = ww.sum()
    if denom == 0:                            # empty after masking -> NaN (matches R 0/0)
        return np.nan                         # empty after masking -> NaN (same as R's 0/0)
    return float((xw * ww).sum() / denom)


def co2_model_reach(
    hrt: float,          # residence time
    Q: float,            # discharge
    C_gw: float,         # groundwater CO2
    C_atm: float,        # atmospheric CO2
    C_bz: float,         # benthic-zone CO2 (calibrated)
    k_co2: float,        # gas-exchange rate
    k_bz: float,         # hyporheic-exchange rate
    F_wc: float,         # water-column reaction (calibrated)
    divergence: int,     # divergent-reach flag
    startflag: int,      # headwater flag
    waterbody: str,      # 'River' | 'Lake/Reservoir' | 'GreatLake'
    upstream_co2: np.ndarray,   # CO2 of directly-upstream reaches
    upstream_Q: np.ndarray,     # Q of directly-upstream reaches
    emergenceQ: float,   # emergent discharge
    C_floor: float | None = None,   # S13: lake floor when F_wc < 0
    benthic_two_way: bool = False,  # v6 option: exchange k_bz (C_bz - C) instead of the one-way source k_bz C_bz
) -> float:
    """Return steady-state CO2 [ppm] for one reach.

    Port of model.R:31-63.
    """
    # --- HEADWATER or DIVERGENT reach: upstream boundary = groundwater --- (model.R:33-36)
    if startflag == 1 or divergence == 2:
        CO2_last = C_gw
        lastQ = 0.0 if emergenceQ >= Q else emergenceQ    # only use emergenceQ if headwater big enough
    else:
        # --- ALL OTHER reaches: discharge-weighted mean of upstream inflow --- (model.R:40-49)
        q_mask = ~np.isnan(upstream_Q)
        if not q_mask.any():
            # no upstream (empty) or all-NA Q -> treat as emergent stream (model.R:44-49)
            CO2_last = C_gw
            lastQ = 0.0 if emergenceQ >= Q else emergenceQ
        else:
            CO2_last = _weighted_mean_narm(upstream_co2, upstream_Q)  # inflow CO2, Q-weighted
            lastQ = float(np.nansum(upstream_Q))                     # total inflow (mass conservation)

    # --- Losing streams: cap inflow Q at this reach's Q so the term cannot exceed it --- (model.R:53)
    if Q < lastQ:
        lastQ = Q

    # --- Great Lakes assumed at atmospheric equilibrium --- (model.R:56-58)
    if waterbody in ("GreatLake", "NonTarget"):   # NonTarget: polygon outside the NLA target population (v6)
        return C_atm

    # --- Backward-difference steady-state solution (equation 6, model.R:60) ---
    numerator = (
        CO2_last
        + ((Q - lastQ) / Q) * C_gw     # lateral groundwater input
        + hrt * k_co2 * C_atm          # gas exchange toward atmosphere
        + hrt * k_bz * C_bz            # hyporheic exchange toward streambed
        + hrt * F_wc                   # net water-column reaction
    )
    denominator = 1.0 + hrt * k_co2 + (Q - lastQ) / Q
    if benthic_two_way:
        denominator = denominator + hrt * k_bz   # two-way form: the benthic term relaxes C toward C_bz
    C = numerator / denominator
    # S13: a negative lake F_wc may not pull the reach below C_floor; raising F_wc to the binding value gives
    # exactly C_floor, so the clip is the per-reach floor rule with the actual inflow.
    if C_floor is not None and F_wc < 0 and waterbody != "River" and C < C_floor:
        C = C_floor
    return C
