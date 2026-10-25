"""
Regression + physics tests for the Python port.

Two kinds of guards
  1. Hand-checkable physics: closed-form / textbook values for the helper functions.
  2. Golden regression: frozen synthetic-basin outputs so any future refactor that
     changes a number is caught. (These are NOT R fixtures, they lock the CURRENT
     Python behavior. Replace/augment with real R output when NAS data is available.)
"""
import numpy as np
import pytest

from conus_carbon import synthetic_basin, calibrate, run_model, hydraulics


# ----------------------------- physics sanity -----------------------------

def test_henry_decreases_with_temperature():
    # Henry's constant falls as water warms (gas less soluble)
    assert float(hydraulics.henry_func(5)) > float(hydraulics.henry_func(25))


def test_henry_positive():
    for T in (0, 10, 20, 30):
        assert float(hydraulics.henry_func(T)) > 0


def test_kco2_steeper_reach_degasses_faster():
    steep = hydraulics.kco2_func(np.array([1.0]), np.array([0.03]), np.array([0.5]),
                                 np.array([15.0]), np.array([0.0]), np.array(["River"]))
    gentle = hydraulics.kco2_func(np.array([0.2]), np.array([0.001]), np.array([1.0]),
                                  np.array([15.0]), np.array([0.0]), np.array(["River"]))
    assert float(steep[0]) > float(gentle[0])


def test_forward_bounded_when_no_internal_sources():
    # With C_bz=0 and F_wc=0, every reach CO2 must lie in [C_atm, C_gw].
    basin = synthetic_basin.make_synthetic_basin()
    h = calibrate.run_forward(basin, (0.0, 0.0, 0.0, 0.0))
    co2 = h["CO2_ppm"].to_numpy()
    assert (co2 >= 400.0 - 1e-6).all() and (co2 <= 16000.0 + 1e-6).all()


def test_cbz_riv_raises_river_median():
    # Identifiable knob: raising Cbz_riv must raise the river median.
    basin = synthetic_basin.make_synthetic_basin()
    t = synthetic_basin.SYNTHETIC_TARGETS
    lo = calibrate.evaluate(basin, (200, 600, 1e-4, 1e-5), **t).river_median
    hi = calibrate.evaluate(basin, (1800, 600, 1e-4, 1e-5), **t).river_median
    assert hi > lo


def test_cross_basin_import_requires_huc4():
    # R always filters imports to the current basin; we require huc4 when downstreamBasin is present.
    import pandas as pd
    basin = synthetic_basin.make_synthetic_basin()
    udf = pd.DataFrame({"downstreamBasin": ["9999"], "exported_ToNode": [int(basin.FromNode.iloc[5])],
                        "exported_CO2_ppm": [9000.0], "exported_Q_cms": [50.0]})
    with pytest.raises(ValueError):
        calibrate.run_forward(basin, (600, 600, 1e-4, 1e-5), upstream_df=udf)  # huc4 missing -> raise


# ----------------------------- golden regression --------------------------
# Frozen CURRENT Python outputs (full precision). If a refactor changes any number,
# these fail and force a conscious decision.

GOLDEN = {
    # params (Cbz_riv, Cbz_lake, Fwc_riv, Fwc_lake): (river_median, lake_median, fitness)
    (600.0, 600.0, 1e-4, 1e-5): (8523.253741755316, 35292.2141379002, 2.456069037339184e-05),
    (0.0, 0.0, 0.0, 0.0):       (6628.180566135038, 794.0778398184458, 0.00022055080362336178),
    (1000.0, 1000.0, 5e-3, 5e-3): (9247.693538900474, 60667.48471716535, 1.4966659164891395e-05),
}


@pytest.mark.parametrize("params,expected", list(GOLDEN.items()))
def test_golden_regression(params, expected):
    basin = synthetic_basin.make_synthetic_basin()
    t = synthetic_basin.SYNTHETIC_TARGETS
    ev = calibrate.evaluate(basin, params, **t)
    assert ev.river_median == pytest.approx(expected[0], rel=0, abs=1e-6)
    assert ev.lake_median == pytest.approx(expected[1], rel=0, abs=1e-6)
    assert ev.fitness == pytest.approx(expected[2], rel=1e-9)
