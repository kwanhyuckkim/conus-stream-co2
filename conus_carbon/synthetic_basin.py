"""
Synthetic test basin generator.

WHY THIS EXISTS
  The real per-basin hydrography comes from NHD-HR .gdb files that live on the
  UMass NAS and are NOT in this repo. To demonstrate the full code flow and to
  prototype the sweep before real data is wired in, we build a small, physically
  plausible toy basin with the same 25-column schema.

THIS IS NOT REAL DATA
  Numbers here are illustrative. The AHG coefficients (a,b,c,f) are placeholder
  literature-typical values because the real ones (depAHG.rds / widAHG.rds) are
  on NAS. Do not treat outputs as scientific results.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import hydraulics

# placeholder at-a-station hydraulic-geometry coefficients (real ones are on NAS)
AHG_A, AHG_B = 2.0, 0.5      # width  W = exp(a)*Q^b  ~ 7.4 * Q^0.5
AHG_C, AHG_F = -1.0, 0.4     # depth  D = exp(c)*Q^f  ~ 0.37 * Q^0.4


def make_synthetic_basin() -> pd.DataFrame:
    """Return a fixed 12-reach toy basin as a 24-column hydrography DataFrame.

    Topology (headwaters at top, single outlet at bottom):

        r1  r2      r3  r4          (headwaters
          \\  \\    /  /
           r5(lake) r6              (mid reaches
             \\    /
              r7(lake)
                |
               r8  r9  r10          more tributaries merging
                \\  |  /
                 r11(lake)
                   |
                  r12               outlet
    """
    # node ids (integers, like NHD FromNode/ToNode)
    # each reach: (reach_id, from_node, to_node, hydroseq, startflag, waterbody, Q, slope, tempC, lenKm, lakeArea_km2)
    rows = [
        # id  from  to   hydroseq start wb                Q     slope  tempC lenKm lakeA
        (1,  101, 201, 120, 1, "River",          0.5,  0.020, 12.0, 1.2, 0.0),
        (2,  102, 201, 119, 1, "River",          0.6,  0.018, 12.5, 1.4, 0.0),
        (3,  103, 202, 118, 1, "River",          0.4,  0.025, 11.0, 1.0, 0.0),
        (4,  104, 202, 117, 1, "River",          0.7,  0.015, 13.0, 1.6, 0.0),
        (5,  201, 203, 100, 0, "Lake/Reservoir", 1.2,  0.001, 12.2, 0.8, 0.30),
        (6,  202, 203,  99, 0, "River",          1.2,  0.012, 12.0, 1.1, 0.0),
        (7,  203, 204,  80, 0, "Lake/Reservoir", 2.6,  0.001, 12.1, 1.5, 1.20),
        (8,  105, 205,  70, 1, "River",          0.8,  0.016, 12.8, 1.3, 0.0),
        (9,  106, 205,  69, 1, "River",          0.5,  0.019, 12.6, 1.0, 0.0),
        (10, 204, 205,  60, 0, "River",          2.7,  0.008, 12.1, 1.7, 0.0),
        (11, 205, 206,  40, 0, "Lake/Reservoir", 4.2,  0.001, 12.3, 2.0, 2.50),
        (12, 206, 207,  20, 0, "River",          4.4,  0.006, 12.2, 2.4, 0.0),
    ]
    cols = ["reach", "FromNode", "ToNode", "HydroSeq", "StartFlag", "waterbody",
            "Q_cms", "Slope", "temp_c", "LengthKM", "lakeArea_km2"]
    df = pd.DataFrame(rows, columns=cols)

    # NHD-like bookkeeping columns
    df["NHDPlusID"] = (900000 + df["reach"]).astype("int64")
    df["WBArea_Permanent_Identifier"] = np.where(
        df["waterbody"] == "Lake/Reservoir", "LK" + df["reach"].astype(str), None
    )
    df["conus"] = 1
    df["Divergence"] = 0
    # StreamOrde is a placeholder (=1) and is NOT used by the model physics/fitness/emissions.
    # A real network would have order rising downstream.
    df["StreamOrde"] = 1
    df["CatchmentAreaSqKm"] = df["Q_cms"] * 50.0     # arbitrary plausible

    # lake volume from area (Cael et al. 2016, setup_hydrography.R:192)
    lake_area_m2 = df["lakeArea_km2"].to_numpy() * 1e6
    lake_vol_m3 = np.where(lake_area_m2 > 0, 0.533 * lake_area_m2 ** 1.204, 0.0)
    # single-throughflow assumption for the toy basin (lakePercent=1)
    df["frac_lakeSurfaceArea_m2"] = lake_area_m2
    df["frac_lakeVol_m3"] = lake_vol_m3

    # AHG coefficients (placeholder)
    a = np.full(len(df), AHG_A); b = np.full(len(df), AHG_B)
    c = np.full(len(df), AHG_C); f = np.full(len(df), AHG_F)

    wb = df["waterbody"].to_numpy()
    Q = df["Q_cms"].to_numpy()
    # derived hydraulics via the ported utils.R functions
    df["D"] = hydraulics.depth_func(wb, Q, df["frac_lakeVol_m3"].to_numpy(),
                                    df["frac_lakeSurfaceArea_m2"].to_numpy(), c, f)
    df["W"] = hydraulics.width_func(wb, Q, a, b)
    df["V"] = hydraulics.velocity_func(wb, Q, df["W"].to_numpy(), df["D"].to_numpy())
    df["HRT"] = hydraulics.restime_water(df["frac_lakeVol_m3"].to_numpy(), df["LengthKM"].to_numpy(),
                                         df["V"].to_numpy(), Q, wb)
    df["henry"] = hydraulics.henry_func(df["temp_c"].to_numpy())
    # NOTE
    # whose area classes are in km2 (0.1/1/10). We replicate that EXACTLY (do not change the
    # flow): any real lake area in m2 exceeds 10, so lakes always take the k600=1.90 class.
    df["k_co2"] = hydraulics.kco2_func(df["V"].to_numpy(), df["Slope"].to_numpy(), df["D"].to_numpy(),
                                       df["temp_c"].to_numpy(), df["frac_lakeSurfaceArea_m2"].to_numpy(), wb)
    df["k_bz"] = hydraulics.kbz_func(df["Slope"].to_numpy(), df["D"].to_numpy(), wb, df["temp_c"].to_numpy())

    from .constants import HYDROGRAPHY_COLUMNS
    return df[list(HYDROGRAPHY_COLUMNS)].copy()


# illustrative GLORICH-style calibration targets for the toy basin (ppm)
SYNTHETIC_TARGETS = {"riverCO2": 2200.0, "lakeCO2": 900.0}
