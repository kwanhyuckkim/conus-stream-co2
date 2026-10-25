"""
Constants for the CONUS CO2 transport model.

All values are taken verbatim from the original R pipeline (`_targets.R`).

Source lines
  _targets.R:38  C_atmosphere <- 400
  _targets.R:39  emergenceQ   <- 0.000105
  _targets.R:40  C_groundwater<- 16000
  _targets.R:43-51 parameter bounds
  _targets.R:54-59 GA meta-parameters
  calibrate_model.R:42 great-lakes skip list
  calibrate_model.R:54-55 photosynthesis HUC2 override
"""
from __future__ import annotations

# ---- Physical constants (fixed, NOT calibrated)
C_ATMOSPHERE_PPM = 400.0        # atmospheric CO2
C_GROUNDWATER_PPM = 16000.0     # groundwater CO2
EMERGENCE_Q_M3S = 0.000105      # emergent headwater discharge

# ---- Calibration parameter bounds
# order
# REFERENCE-ONLY: the R GA uses these to bound its search. This port has no calibration/sweep
# driver yet, so nothing here is enforced. A future sweep MUST apply these bounds (and the
# photosynthesis override below) when it searches parameters.
LOWER_BOUNDS = (0.0, 0.0, 0.0, 0.0)                 # [ppm, ppm, ppm/s, ppm/s]
UPPER_BOUNDS = (C_GROUNDWATER_PPM, C_GROUNDWATER_PPM, 1e-2, 1e-2)  # [ppm, ppm, ppm/s, ppm/s]
# Cbz ceilings track the groundwater endmember (bed pore water is a stream-groundwater
# mixture, so it cannot exceed C_gw); replaces the undocumented 2,000 constant.
PARAM_NAMES = ("Cbz_riv", "Cbz_lake", "Fwc_riv", "Fwc_lake")

# Great Lakes / Great Basin HUC2 prefixes allow net lake photosynthesis (negative Fwc_lake). (calibrate_model.R:54-55)
# NOT wired into any driver yet, a sweep on HUC2 04/16 must widen Fwc_lake to [lower, upper] below.
PHOTOSYNTHESIS_HUC2 = ("04", "16")
FWC_LAKE_LOWER_PHOTO = -1e-5
FWC_LAKE_UPPER_PHOTO = 0.0

# ---- Special basin lists
GREAT_LAKES_HUC4 = ("0418", "0419", "0424", "0426", "0428")   # calibrate_model.R:42 / run_model.R:26
INDIANA_HUC4 = ("0508", "0509", "0514", "0512", "0712", "0404", "0405", "0410")  # setup_hydrography.R:21

# ---- GA meta-parameters (reference only; we are replacing the GA) ---- (_targets.R:54-59)
GA_POP_SIZE = 25
GA_MAX_ITER = 500
GA_RUN = 50            # early-stop after 50 stagnant generations
GA_MAX_FITNESS = 1.0 / 10.0   # [1/ppm]
GA_MUTATION_RATE = 0.10
GA_SEED = 12

# ---- Output hydrography schema (24 unique columns)
# NOTE: R's select() lists 'NHDPlusID' twice, so it reads as "25" but is 24 unique columns.
HYDROGRAPHY_COLUMNS = (
    "NHDPlusID", "WBArea_Permanent_Identifier", "conus", "FromNode", "ToNode",
    "HydroSeq", "StartFlag", "Divergence", "waterbody", "StreamOrde", "LengthKM",
    "Slope", "Q_cms", "frac_lakeSurfaceArea_m2", "frac_lakeVol_m3",
    "CatchmentAreaSqKm", "temp_c", "HRT", "henry", "k_co2", "k_bz", "W", "D", "V",
)
