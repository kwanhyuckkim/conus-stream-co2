"""
Build the per-basin hydrography routing table from NHD-HR.

STRUCTURAL PORT
  `src/setup_hydrography.R` reads NHD-HR .gdb geodatabases (NHDFlowline,
  NHDWaterbody, NHDPlusEROMMA, NHDPlusFlowlineVAA, NHDPlusIncrTempMM01-12) plus
  hydraulic-geometry models (depAHG.rds / widAHG.rds). Those files live on the
  UMass NAS and are NOT in this repo, so this module documents the exact pipeline
  and provides the column contract, but the .gdb reader is left as a TODO
  (use pyogrio when the data is mounted). For now, use `synthetic_basin` to
  exercise the flow.

PIPELINE STEPS (setup_hydrography.R)
  1. Read NHDFlowline (+ Indiana shp override) and NHDWaterbody. (R:28-47)
  2. Filter lakes to FType in {390, 436}. (R:50-52)
  3. Read EROM (mean annual flow), VAA, and 12 monthly temperature layers. (R:53-137)
  4. Mean annual air temp = row mean of monthly temps. (R:139-141)
  5. Join lakes/EROM/VAA/temp onto flowlines; fix duplicated column names. (R:143-161)
  6. Flag CONUS reaches by intersecting the US boundary shapefile. (R:163-167)
  7. StreamOrde = StreamCalc; Q_cms = QBMA * 0.0283 (cfs -> cms). (R:170-171)
  8. waterbody = 'Lake/Reservoir' if it has a lake area, else 'River'. (R:181)
  9. Filter StreamOrde>0, HydroSeq not NA, FlowDir==1, Q_cms>0. (R:189-190)
 10. lakeVol_m3 = 0.533 * (LakeAreaSqKm*1e6)^1.204  (Cael 2016). (R:192)
 11. Fix bad (<=0) slopes and missing temps from neighbors. (R:194-208)
 12. Fractional lake area/volume per throughflow line (Brinkerhoff 2021). (R:213-220)
 13. AHG coefficients a,b (width) and c,f (depth) from *.rds. (R:224-229)
 14. D=depth_func, W=width_func, V=velocity_func. (R:232-234)
 15. temp_c = 3.941 + 0.818*airTemp (Lauerwald 2015 air->water). (R:237)
 16. HRT, henry, k_co2, k_bz via the ported hydraulics functions. (R:240-247)
 17. Select the 25-column output schema. (R:255-256)
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import hydraulics
from .constants import HYDROGRAPHY_COLUMNS, INDIANA_HUC4

# Cael et al. 2016 lake volume-from-area coefficients
CAEL_A, CAEL_B = 0.533, 1.204
# Lauerwald et al. 2015 air->water temperature regression
LAUERWALD_INTERCEPT, LAUERWALD_SLOPE = 3.941, 0.818


def air_to_water_temp(air_temp_c):
    """Water temp from mean annual air temp (Lauerwald 2015). setup_hydrography.R:237."""
    return LAUERWALD_INTERCEPT + LAUERWALD_SLOPE * np.asarray(air_temp_c)


def lake_volume_m3(lake_area_sqkm):
    """Lake volume [m3] from surface area [km2] (Cael 2016). setup_hydrography.R:192."""
    return CAEL_A * (np.asarray(lake_area_sqkm) * 1e6) ** CAEL_B


def derive_model_columns(nhd: pd.DataFrame, ahg: dict) -> pd.DataFrame:
    """Given a joined/cleaned flowline table + AHG coefficients, compute the model columns.

    `ahg` keys: a, b (width), c, f (depth). `nhd` must already contain the raw
    attributes: waterbody, Q_cms, Slope, LengthKM, frac_lakeVol_m3,
    frac_lakeSurfaceArea_m2, temp_c (water temp).
    """
    nhd = nhd.copy()
    wb = nhd["waterbody"].to_numpy()
    Q = nhd["Q_cms"].to_numpy(dtype=float)
    a, b, c, f = ahg["a"], ahg["b"], ahg["c"], ahg["f"]

    nhd["D"] = hydraulics.depth_func(wb, Q, nhd["frac_lakeVol_m3"].to_numpy(),
                                     nhd["frac_lakeSurfaceArea_m2"].to_numpy(), c, f)
    nhd["W"] = hydraulics.width_func(wb, Q, a, b)
    nhd["V"] = hydraulics.velocity_func(wb, Q, nhd["W"].to_numpy(), nhd["D"].to_numpy())
    nhd["HRT"] = hydraulics.restime_water(nhd["frac_lakeVol_m3"].to_numpy(),
                                          nhd["LengthKM"].to_numpy(), nhd["V"].to_numpy(), Q, wb)
    nhd["henry"] = hydraulics.henry_func(nhd["temp_c"].to_numpy())
    # lakeArea passed in m2 to match setup_hydrography.R:246 exactly
    nhd["k_co2"] = hydraulics.kco2_func(nhd["V"].to_numpy(), nhd["Slope"].to_numpy(),
                                        nhd["D"].to_numpy(), nhd["temp_c"].to_numpy(),
                                        nhd["frac_lakeSurfaceArea_m2"].to_numpy(), wb)
    nhd["k_bz"] = hydraulics.kbz_func(nhd["Slope"].to_numpy(), nhd["D"].to_numpy(),
                                      wb, nhd["temp_c"].to_numpy())
    return nhd


def setup_hydrography(path_to_data: str, huc4: str) -> pd.DataFrame:
    """Read + build + clean the NHD-HR hydrography into the routing table.

    TODO: implement the .gdb reader with pyogrio once NAS data is mounted.
    """
    raise NotImplementedError(
        "setup_hydrography needs NHD-HR .gdb files (on NAS) + AHG .rds coefficients. "
        "Use conus_carbon.synthetic_basin.make_synthetic_basin() to exercise the flow, "
        "or wire pyogrio to read '{path}/HUC2_{h2}/NHDPLUS_H_{huc4}_HU4_GDB/...'.".format(
            path=path_to_data, h2=huc4[:2], huc4=huc4
        )
    )


__all__ = ["setup_hydrography", "derive_model_columns", "air_to_water_temp",
           "lake_volume_m3", "HYDROGRAPHY_COLUMNS", "INDIANA_HUC4"]
