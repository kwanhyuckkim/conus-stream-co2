"""
conus_carbon: Python port of the CONUS stream-CO2 reactive-transport model.

Faithful port of the R pipeline (Saccardi et al. 2024, AGU Advances). The forward
model, the reach routing, and the calibration fitness reproduce the R code's flow
exactly; only implementation efficiency (numpy vectorization, an O(N) upstream-
index map) differs, with identical numerical results.

Quick start
    from conus_carbon import synthetic_basin, calibrate
    basin = synthetic_basin.make_synthetic_basin()
    ev = calibrate.evaluate(basin, params=(1000, 1000, 5e-3, 5e-3),
                            riverCO2=2200, lakeCO2=900)
    print(ev.river_median, ev.lake_median, ev.fitness)
"""
from __future__ import annotations

from . import calibrate, constants, hydraulics, model, run_model, synthetic_basin
from .calibrate import Evaluation, basin_medians, evaluate, run_forward
from .model import co2_model_reach
from .run_model import run_model as run_calibrated_model

__version__ = "0.1.0"

__all__ = [
    "constants", "hydraulics", "model", "calibrate", "run_model", "synthetic_basin",
    "co2_model_reach", "run_forward", "basin_medians", "evaluate", "Evaluation",
    "run_calibrated_model",
]
