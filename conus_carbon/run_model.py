"""
Run the calibrated model and compute per-reach CO2 flux to the atmosphere.

Faithful port of `src/run_model.R` (runModel).

The reach loop is identical to calibrate.run_forward; here we reuse it and then
attach the diagnostic + emission output columns exactly as run_model.R does.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .calibrate import run_forward
from .constants import C_ATMOSPHERE_PPM, C_GROUNDWATER_PPM, EMERGENCE_Q_M3S, GREAT_LAKES_HUC4
from .hydraulics import schmidt_co2

# output schema of run_model.R (lines 40 / 107)
OUTPUT_COLUMNS = (
    "NHDPlusID", "WBArea_Permanent_Identifier", "ToNode", "FromNode", "conus",
    "CO2_ppm", "StreamOrde", "k_co2_m_s", "k600_m_s", "LengthKM", "Q_m3_s", "W_m",
    "lakeSA_m2", "waterbody", "V_m_s", "Slope", "Water_temp_c", "henry", "FCO2_gC_m2_yr",
)


def _add_output_columns(h: pd.DataFrame, C_atm: float) -> pd.DataFrame:
    """Attach diagnostic + emission columns. Port of run_model.R:29-38 / 96-105."""
    h = h.copy()
    temp = h["temp_c"].to_numpy(dtype=float)
    # Schmidt number to convert k_co2 back to a k600 (diagnostic); single source in hydraulics.
    sc = schmidt_co2(temp)
    h["HRT_s"] = h["HRT"]                                   # [s]
    h["k_co2_m_s"] = h["k_co2"] * h["D"]                    # [m/s] = [1/s] * depth
    h["k600_m_s"] = h["k_co2_m_s"] * (600.0 / sc) ** -0.5   # [m/s] back to k600 (diagnostic)
    h["k_bz_m_s"] = h["k_bz"] * h["D"]                      # [m/s]
    h["W_m"] = h["W"]                                       # [m]
    h["V_m_s"] = h["V"]                                     # [m/s]
    h["Q_m3_s"] = h["Q_cms"]                                # [m3/s]
    h["lakeSA_m2"] = h["frac_lakeSurfaceArea_m2"]           # [m2]
    h["Water_temp_c"] = h["temp_c"]                         # [c]
    # CO2 flux to atmosphere [g-C/m2/yr] (run_model.R:38,105)
    h["FCO2_gC_m2_yr"] = (
        h["k_co2_m_s"]
        * (((h["CO2_ppm"] - C_atm) * h["henry"]) / 1_000_000.0)  # ppm*henry -> mol/L /1e6
        * (1.0 / 0.001)                                          # L -> m3
        * 12.01                                                  # mol C -> g C
        * (60 * 60 * 24 * 365)                                   # per second -> per year
    )
    return h


def run_model(
    hydrography: pd.DataFrame,
    calibration_results: dict,      # keys: Cbz_riv, Cbz_lake, Fwc_riv, Fwc_lake
    C_gw: float = C_GROUNDWATER_PPM,
    C_atm: float = C_ATMOSPHERE_PPM,
    huc4: str = "",
    emergenceQ: float = EMERGENCE_Q_M3S,
    upstream_df: pd.DataFrame | None = None,
    benthic_two_way: bool = False,
) -> pd.DataFrame:
    """Run the calibrated CO2 model for a basin. Port of run_model.R:24-112."""
    # Great Lakes basins: force atmospheric equilibrium (run_model.R:26-27)
    if huc4 in GREAT_LAKES_HUC4:
        h = hydrography.copy()
        h["CO2_ppm"] = C_atm
    else:
        params = (
            calibration_results["Cbz_riv"], calibration_results["Cbz_lake"],
            calibration_results["Fwc_riv"], calibration_results["Fwc_lake"],
        )
        h = run_forward(hydrography, params, C_gw=C_gw, C_atm=C_atm,
                        emergenceQ=emergenceQ, upstream_df=upstream_df, huc4=huc4,
                        benthic_two_way=benthic_two_way)

    h = _add_output_columns(h, C_atm)
    return h[list(OUTPUT_COLUMNS)]
