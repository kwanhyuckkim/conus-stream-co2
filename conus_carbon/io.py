"""
Loading helpers: calibration targets (in the repo) and a prepared basin (from the lab).
"""
from __future__ import annotations

import numbers
from pathlib import Path

import pandas as pd

from .constants import HYDROGRAPHY_COLUMNS

# repo data folder
DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_calibration_targets(huc4: str, path: str | Path | None = None) -> dict:
    """Return {'riverCO2', 'lakeCO2'} GLORICH targets [ppm] for a HUC4 from HUC4_calibration.csv.

    This file IS in the repo, so targets are available now (unlike the hydrography)."""
    path = Path(path) if path is not None else DATA_DIR / "HUC4_calibration.csv"
    df = pd.read_csv(path, dtype={"HUC4": str, "HUC2": str})
    # numbers.Number covers python int/float AND numpy scalars (np.int64/np.float64), mirroring
    # sweep._norm_huc4 so int 418 / np.float64(418.0) both resolve to '0418'.
    key = str(int(huc4)) if isinstance(huc4, numbers.Number) else str(huc4)
    if key.isdigit():
        key = key.zfill(4)                 # mirror sweep._norm_huc4 so int 418 -> '0418'
    row = df[df["HUC4"] == key]
    if row.empty:
        raise KeyError(f"HUC4 {huc4!r} not found in {path}")
    river = float(row["River"].iloc[0])
    lake = float(row["Lake"].iloc[0])
    # a blank/NaN target would silently NaN-poison every cost downstream; fail loudly instead.
    import numpy as np
    if not (np.isfinite(river) and np.isfinite(lake)):
        raise ValueError(f"HUC4 {huc4!r} has a non-finite target (River={river}, Lake={lake}) in {path}")
    return {"riverCO2": river, "lakeCO2": lake}


def load_basin_hydrography(path: str | Path) -> pd.DataFrame:
    """Load a prepared 24-column basin hydrography from CSV or Parquet.

    Use this once the lab shares a real basin (the output of setup_hydrography). It validates
    the 24-column schema so a mislabeled file fails loudly instead of producing wrong numbers."""
    path = Path(path)
    if path.suffix.lower() in (".parquet", ".pq"):
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path)
    missing = [c for c in HYDROGRAPHY_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required hydrography columns: {missing}. "
            f"Expected the 24-column schema from setup_hydrography (conus_carbon.constants.HYDROGRAPHY_COLUMNS)."
        )
    cols = list(HYDROGRAPHY_COLUMNS)
    if "C_gw_reach" in df.columns:
        cols.append("C_gw_reach")
    if "C_atm_reach" in df.columns:
        cols.append("C_atm_reach")
    return df[cols].copy()
