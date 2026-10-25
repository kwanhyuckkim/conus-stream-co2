#!/usr/bin/env python3
"""
Run the Python-ported model on a prepared hydrography exported from R, and save per-reach output.

This isolates exactly what was ported (CO2_MODEL + reach loop + emissions). It does NOT re-run
setup_hydrography (that needs the NHD-HR .gdb and is not ported); instead it consumes the SAME
24-column hydrography that R used, so any difference is purely model logic.

Usage:
  python run_python_model.py --hydrography <huc4>_hydrography.csv \
                             --params <huc4>_params.csv \
                             --huc4 0107 --out out/python/0107_py_output.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# make the conus_carbon package importable from the repo root
_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from conus_carbon import io, run_model                       # noqa: E402
from conus_carbon.constants import C_ATMOSPHERE_PPM, C_GROUNDWATER_PPM, EMERGENCE_Q_M3S  # noqa: E402


def load_params(path: Path) -> tuple[dict, dict]:
    """Read the one-row params CSV. Returns (knobs, constants).

    The R export writes Cgw/Catm/emergenceQ into this file ON PURPOSE, so BOTH models run identical
    physical constants. We read them here (falling back to constants.py only if a column is absent,
    with a warning), so R's export is the single source of truth for the comparison."""
    df = pd.read_csv(path)
    row = df.iloc[0]
    knobs = {k: float(row[k]) for k in ("Cbz_riv", "Cbz_lake", "Fwc_riv", "Fwc_lake")}
    consts = {}
    for col, key, default in (("Cgw", "C_gw", C_GROUNDWATER_PPM),
                              ("Catm", "C_atm", C_ATMOSPHERE_PPM),
                              ("emergenceQ", "emergenceQ", EMERGENCE_Q_M3S)):
        if col in df.columns and pd.notna(row[col]):
            consts[key] = float(row[col])
        else:
            consts[key] = float(default)
            print(f"[python] WARNING: '{col}' not in params.csv; using constants.py default {default}. "
                  f"Comparison then assumes R used the same value.", file=sys.stderr)
    return knobs, consts


def main():
    ap = argparse.ArgumentParser(description="Run the Python-port model on an R-exported hydrography.")
    ap.add_argument("--hydrography", required=True, help="<huc4>_hydrography.csv (24-col schema)")
    ap.add_argument("--params", required=True, help="<huc4>_params.csv (Cbz_riv,Cbz_lake,Fwc_riv,Fwc_lake)")
    ap.add_argument("--huc4", required=True, help="basin id, e.g. 0107")
    ap.add_argument("--out", required=True, help="output CSV path")
    # constants come from params.csv (R's export). CLI flags override only if you know why.
    ap.add_argument("--cgw", type=float, default=None, help="override groundwater CO2 [ppm]")
    ap.add_argument("--catm", type=float, default=None, help="override atmospheric CO2 [ppm]")
    ap.add_argument("--emergenceq", type=float, default=None, help="override emergent Q [m3/s]")
    args = ap.parse_args()

    hydro = io.load_basin_hydrography(args.hydrography)       # validates the 24-column schema
    knobs, consts = load_params(Path(args.params))
    # params.csv drives the constants; a CLI flag overrides only if explicitly given
    C_gw = args.cgw if args.cgw is not None else consts["C_gw"]
    C_atm = args.catm if args.catm is not None else consts["C_atm"]
    emergenceQ = args.emergenceq if args.emergenceq is not None else consts["emergenceQ"]

    out = run_model.run_model(
        hydro, knobs, C_gw=C_gw, C_atm=C_atm, huc4=str(args.huc4),
        emergenceQ=emergenceQ, upstream_df=None,               # single-basin validation (no imports)
    )
    # keep the columns the comparison needs (reach id + the two modeled quantities)
    keep = out[["NHDPlusID", "CO2_ppm", "FCO2_gC_m2_yr"]].copy()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    keep.to_csv(args.out, index=False)
    n_ok = int(keep["CO2_ppm"].notna().sum())
    print(f"[python] HUC4 {args.huc4}: {len(keep)} reaches ({n_ok} finite CO2) -> {args.out}")


if __name__ == "__main__":
    main()
