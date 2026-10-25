"""
Hydraulic + gas-exchange helper functions.

Faithful port of the physics helpers in `src/utils.R`.

Ported functions
  henry_func      (15)  Henry's constant
  restime_water   (43)  residence time
  width_func      (74)  river width (AHG)
  depth_func      (100) depth (AHG
  velocity_func   (124) river velocity
  kco2_func       (152) gas-exchange rate
  kbz_func        (192) hyporheic-exchange rate

EFFICIENCY
mapply). Same math, same numbers, applied to a whole column at once.
"""
from __future__ import annotations

import numpy as np


def schmidt_co2(temp):
    """Schmidt number for CO2 (Raymond 2012 / Wanninkhof 1991). Same formula used by kco2 and kbz.
    Public so run_model can reuse the single source of these coefficients."""
    return 1911.0 - 118.11 * temp + 3.453 * temp ** 2 - 0.0413 * temp ** 3


def henry_func(temp):
    """Temperature-dependent Henry's law constant [mol/(L*atm)]. utils.R:15-26."""
    A, B, C, D, E = 108.3865, 0.01985076, -6919.53, -40.4515, 669365.0   # constants
    temp_k = np.asarray(temp) + 273.15                                    # to Kelvin
    return 10.0 ** (A + B * temp_k + C / temp_k + D * np.log10(temp_k) + E / temp_k ** 2)


def restime_water(Vol, lengthKm, v, Q, waterbody):
    """Water residence time [s]. Lake: Vol/Q. River: length/velocity. utils.R:43-58."""
    Vol, lengthKm, v, Q = map(np.asarray, (Vol, lengthKm, v, Q))
    waterbody = np.asarray(waterbody)
    out = np.full(np.broadcast(Vol, lengthKm, v, Q).shape, np.nan, dtype=float)
    nonzero = Q != 0                                    # Q==0 -> NA (R:44-45)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(nonzero & (waterbody == "Lake/Reservoir"), Vol / Q, out)          # lake: Vol/Q
        out = np.where(nonzero & (waterbody == "River"), (lengthKm * 1000.0) / v, out)   # river: length_m/velocity
    return out


def width_func(waterbody, Q, a, b):
    """River width via at-a-station hydraulic geometry W = exp(a)*Q^b [m]. Lakes -> NaN. utils.R:74-82."""
    Q = np.asarray(Q)
    with np.errstate(invalid="ignore"):
        return np.where(np.asarray(waterbody) == "River", np.exp(a) * Q ** b, np.nan)


def depth_func(waterbody, Q, lakeVol, lakeArea, c, f):
    """Depth [m]. River: exp(c)*Q^f. Lake: volume/area. utils.R:100-108."""
    Q, lakeVol, lakeArea = map(np.asarray, (Q, lakeVol, lakeArea))
    with np.errstate(divide="ignore", invalid="ignore"):
        river_d = np.exp(c) * Q ** f
        lake_d = lakeVol / lakeArea
    return np.where(np.asarray(waterbody) == "River", river_d, lake_d)


def velocity_func(waterbody, Q, W, D):
    """River velocity via continuity V = Q/(W*D) [m/s]. Lakes -> NaN. utils.R:124-132."""
    Q, W, D = map(np.asarray, (Q, W, D))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.asarray(waterbody) == "River", Q / (W * D), np.nan)


def kco2_func(velocity, slope, depth, temp, lakeArea, waterbody):
    """CO2 gas-exchange rate constant k_co2 [1/s]. utils.R:152-175.

    River: Ulseth (2019) two-regime k600 from energy dissipation eD = 9.8*V*S.
    Lake : Read (2012)/Raymond (2013) k600 by surface-area class.

    UNIT QUIRK (kept faithful)
      The area classes below are in km2, but setup_hydrography.R:246 passes lakeArea
      in m2. We do NOT change this. Consequence: any real lake (m2 >> 10) always
      lands in the last class (k600 = 1.90 m/day); a very small lake split across many
      throughflow lines can have a fractional area < 10 m2 and fall in a lower class.
    """
    velocity, slope, depth, temp, lakeArea = map(
        np.asarray, (velocity, slope, depth, temp, lakeArea)
    )
    waterbody = np.asarray(waterbody)
    sc = schmidt_co2(temp)                                    # Schmidt number

    with np.errstate(divide="ignore", invalid="ignore"):
        eD = 9.8 * velocity * slope                           # energy dissipation [m2/s3]
        # two-regime river k600 (Ulseth 2019)
        k600_river = np.where(eD <= 0.02,
                              np.exp(3.10 + 0.35 * np.log(eD)),
                              np.exp(6.43 + 1.18 * np.log(eD)))   # [m/day]
        # lake k600 by area class [km2] (utils.R:163-166)
        k600_lake = np.where(lakeArea <= 0.1, 0.54,
                     np.where(lakeArea <= 1.0, 1.16,
                       np.where(lakeArea <= 10.0, 1.32, 1.90)))    # [m/day]
        k600 = np.where(waterbody == "River", k600_river, k600_lake)
        # convert k600 -> kco2 and to [1/s] (utils.R:171-172)
        k_co2 = k600 / (600.0 / sc) ** -0.5                   # Schmidt correction
        k_co2 = (k_co2 / depth) / (24 * 60 * 60)              # per depth, per second
    return k_co2


def kbz_func(slope, depth, waterbody, temp):
    """Hyporheic/benthic exchange rate constant k_bz [1/s]. utils.R:192-206.

    River: Grant (2018). Lake: Lorke & Peeters (2006). Both scale with shear velocity.
    """
    slope, depth, temp = map(np.asarray, (slope, depth, temp))
    waterbody = np.asarray(waterbody)
    sc = schmidt_co2(temp)                                   # Schmidt number
    with np.errstate(invalid="ignore"):
        Ustar = np.sqrt(9.8 * depth * slope)                 # shear velocity [m/s]
        kbz_river = 0.3 * Ustar * sc ** (-2.0 / 3.0)         # Grant 2018 [m/s]
        kbz_lake = (1.0 / 9.0) * sc ** (-1.0 / 2.0) * Ustar  # Lorke & Peeters 2006 [m/s]
        kbz = np.where(waterbody == "River", kbz_river, kbz_lake)
        kbz = kbz / depth                                    # to [1/s]
    return kbz
