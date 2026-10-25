"""
Forward simulation over a basin + the calibration fitness.

Faithful port of `src/calibrate_model.R` (calibrateModel, lines 127-197).

FLOW IS UNCHANGED
  reaches are sorted upstream->downstream (descending HydroSeq) and a sequential
  loop fills each reach's CO2, exactly like the R code. The fitness is the same:
  median modeled pCO2 per water-body vs the two GLORICH targets, cost summed,
  fitness = 1/cost.

EFFICIENCY-ONLY CHANGE
  R searched `which(toNode_vec == fromNode)` inside every reach (O(N^2) total).
  We precompute an upstream-index map once (O(N)); results are identical.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .constants import C_ATMOSPHERE_PPM, C_GROUNDWATER_PPM, EMERGENCE_Q_M3S
from .model import co2_model_reach


def build_upstream_map(from_node: np.ndarray, to_node_full: np.ndarray):
    """upstream_map[i] = indices j where to_node[j] == from_node[i] (reaches flowing INTO reach i).

    This replaces the per-reach `which()` scan with one O(N) precomputation.
    """
    # NaN-safe: R's NA node never matches (which(NA == fromNode) drops it). Skip NaN instead of
    # casting to int64, which would emit a RuntimeWarning + a platform-defined sentinel value.
    idx_by_tonode: dict[int, list[int]] = defaultdict(list)
    for j, tn in enumerate(to_node_full):
        if np.isfinite(tn):
            idx_by_tonode[int(tn)].append(j)
    out = []
    for fn in from_node:
        out.append(np.asarray(idx_by_tonode.get(int(fn), []), dtype=int)
                   if np.isfinite(fn) else np.asarray([], dtype=int))
    return out


def run_forward(
    hydrography: pd.DataFrame,
    params,                       # (Cbz_riv, Cbz_lake, Fwc_riv, Fwc_lake)
    C_gw: float = C_GROUNDWATER_PPM,
    C_atm: float = C_ATMOSPHERE_PPM,
    emergenceQ: float = EMERGENCE_Q_M3S,
    upstream_df: pd.DataFrame | None = None,
    huc4: str | None = None,      # needed to filter cross-basin imports
    benthic_two_way: bool = False,  # v6 option, see model.co2_model_reach
) -> pd.DataFrame:
    """Run the reach loop and return the sorted hydrography with a new `CO2_ppm` column.
    """
    Cbz_riv, Cbz_lake, Fwc_riv, Fwc_lake = params

    # sort reaches upstream -> downstream (HydroSeq descending)
    # R's filter(HydroSeq != 0) drops NA HydroSeq too (NA -> FALSE); replicate with notna() & != 0.
    # HydroSeq is a unique per-reach routing key, so tie order never arises and results are
    # order-independent (reaches sharing a HydroSeq would be parallel, not upstream of each other).
    # kind="stable" is used only to make the sort deterministic; it does NOT byte-match R's
    # order(-HydroSeq) tie order (pandas reverses ties under ascending=False), which is immaterial here.
    h = hydrography[hydrography["HydroSeq"].notna() & (hydrography["HydroSeq"] != 0)].sort_values(
        "HydroSeq", ascending=False, kind="stable"
    ).reset_index(drop=True)
    n = len(h)

    # per-reach fixed inputs (precomputed in setup)
    start = h["StartFlag"].to_numpy()
    diverg = h["Divergence"].to_numpy()
    from_node = h["FromNode"].to_numpy()
    to_node = h["ToNode"].to_numpy()
    waterbody = h["waterbody"].to_numpy()
    hrt = h["HRT"].to_numpy(dtype=float)
    kco2 = h["k_co2"].to_numpy(dtype=float)
    kbz = h["k_bz"].to_numpy(dtype=float)
    Q = h["Q_cms"].to_numpy(dtype=float)

    CO2 = np.full(n, np.nan)
    cgw_vec = (h["C_gw_reach"].to_numpy(dtype=float)
               if "C_gw_reach" in h.columns else np.full(n, C_gw))
    catm_vec = (h["C_atm_reach"].to_numpy(dtype=float)
                if "C_atm_reach" in h.columns else np.full(n, C_atm))   # lake effective equilibrium (S11)                     # results vector

    # S13: when the lake F_wc knob is negative, the reach solver raises F_wc per reach so that no lake reach
    # falls below C_FLOOR (NLA 2012 5th percentile, 57 uatm), given its actual inflow.
    C_FLOOR = 57.0 if Fwc_lake < 0 else None
    # map the 4 calibrated knobs to per-reach values
    bz = np.where(waterbody == "River", Cbz_riv, Cbz_lake)     # C_bz: river vs lake
    Fwc = np.where(waterbody == "River", Fwc_riv, Fwc_lake)    # F_wc: river vs lake

    # append cross-basin upstream imports to the search arrays
    if upstream_df is not None and len(upstream_df) > 0:
        # R ALWAYS filters the export table to THIS basin (calibrate_model.R:161 / run_model.R:70),
        # so we require huc4 when a downstreamBasin column is present rather than silently skipping it.
        if "downstreamBasin" in upstream_df.columns:
            if huc4 is None:
                raise ValueError(
                    "upstream_df has a 'downstreamBasin' column but huc4 was not provided; "
                    "R always filters imports to the current basin (pass huc4=...)."
                )
            upstream_df = upstream_df[upstream_df["downstreamBasin"] == huc4]
    if upstream_df is not None and len(upstream_df) > 0:
        to_node_full = np.concatenate([to_node, upstream_df["exported_ToNode"].to_numpy()])
        CO2_full = np.concatenate([CO2, upstream_df["exported_CO2_ppm"].to_numpy()])
        Q_full = np.concatenate([Q, upstream_df["exported_Q_cms"].to_numpy()])
    else:
        to_node_full, CO2_full, Q_full = to_node, CO2, Q

    # precompute upstream index map (the O(N) optimization)
    upstream_map = build_upstream_map(from_node, to_node_full)

    # sequential reach loop, upstream -> downstream
    for i in range(n):
        up = upstream_map[i]
        CO2_full[i] = co2_model_reach(
            hrt=hrt[i], Q=Q[i], C_gw=cgw_vec[i], C_atm=catm_vec[i],
            C_bz=bz[i], k_co2=kco2[i], k_bz=kbz[i], F_wc=Fwc[i],
            divergence=diverg[i], startflag=start[i], waterbody=waterbody[i],
            upstream_co2=CO2_full[up], upstream_Q=Q_full[up], emergenceQ=emergenceQ, C_floor=C_FLOOR,
            benthic_two_way=benthic_two_way,
        )

    co2_basin = CO2_full[:n]                      # drop appended rows
    h = h.copy()
    h["CO2_ppm"] = co2_basin
    return h


def _grouped_medians(h: pd.DataFrame):
    """Median CO2_ppm per water-body (skips NaN, like R median(na.rm=T)). Shared by the two callers."""
    return h.groupby("waterbody")["CO2_ppm"].median()


def basin_medians(hydrography_with_co2: pd.DataFrame) -> tuple[float, float]:
    """Return (median river pCO2, median lake pCO2) over the basin."""
    med = _grouped_medians(hydrography_with_co2)
    river = float(med["River"]) if "River" in med.index else np.nan
    lake = float(med["Lake/Reservoir"]) if "Lake/Reservoir" in med.index else np.nan
    return river, lake


def cost_from_medians(river_median, lake_median, riverCO2, lakeCO2,
                      has_river=True, has_lake=True) -> float:
    """cost = |median_lake - lakeCO2| + |median_river - riverCO2|; absent water-body drops its term."""
    terms = []
    if has_lake:
        terms.append(abs(lake_median - lakeCO2))
    if has_river:
        terms.append(abs(river_median - riverCO2))
    return float(np.sum(terms))


@dataclass
class Evaluation:
    """One fitness evaluation's outputs (everything the sweep needs)."""
    river_median: float
    lake_median: float
    cost: float
    fitness: float          # 1/cost (R converts to maximization)


def evaluate(
    hydrography: pd.DataFrame,
    params,
    riverCO2: float,
    lakeCO2: float,
    huc4: str | None = None,      # forwarded so cross-basin imports are filtered like R
    **kwargs,
) -> Evaluation:
    """Run forward + compute medians, cost, and fitness for one parameter vector.
    """
    h = run_forward(hydrography, params, huc4=huc4, **kwargs)
    med = _grouped_medians(h)
    has_river = "River" in med.index
    has_lake = "Lake/Reservoir" in med.index
    river_med = float(med["River"]) if has_river else np.nan
    lake_med = float(med["Lake/Reservoir"]) if has_lake else np.nan
    cost = cost_from_medians(river_med, lake_med, riverCO2, lakeCO2, has_river, has_lake)
    fitness = 1.0 / cost if cost != 0 else np.inf     # R: out <- 1/cost
    return Evaluation(river_median=river_med, lake_median=lake_med, cost=cost, fitness=fitness)
