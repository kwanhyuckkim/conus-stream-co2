#!/usr/bin/env python3
"""
Compare the R reference output against the Python-port output, reach by reach.

Inputs (per basin)
  <huc4>_R_output.csv   : from validation/R/export_reference.R (columns: NHDPlusID, CO2_ppm, FCO2_gC_m2_yr)
  <huc4>_py_output.csv  : from validation/python/run_python_model.py (same columns)

It merges on NHDPlusID and reports, per basin and in aggregate:
  - n reaches, n matched
  - max / median absolute error on CO2_ppm and FCO2
  - max relative error
  - number of reaches exceeding the tolerance
  - PASS / FAIL at the tolerance

Exit code 0 if every basin passes, 1 otherwise (usable in a SLURM chain).

Usage:
  python compare.py --ref-dir out/R --py-dir out/python --huc4 0107 0108 ... \
                    --abs-tol 1e-6 --rel-tol 1e-6 --report out/compare_report.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KEY = "NHDPlusID"
COMPARE_COLS = ["CO2_ppm", "FCO2_gC_m2_yr"]   # the quantities that must match R


def _both_nan(a, b):
    return np.isnan(a) & np.isnan(b)


def compare_basin(ref_csv: Path, py_csv: Path, abs_tol: float, rel_tol: float) -> dict:
    """Compare one basin's R vs Python output. Returns a summary dict."""
    ref = pd.read_csv(ref_csv)
    py = pd.read_csv(py_csv)
    for name, df in (("R", ref), ("Python", py)):
        missing = [c for c in [KEY] + COMPARE_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"{name} output {ref_csv if name=='R' else py_csv} missing columns: {missing}")

    # merge on reach ID; 'outer' so we also catch reaches present in one but not the other
    m = ref[[KEY] + COMPARE_COLS].merge(
        py[[KEY] + COMPARE_COLS], on=KEY, how="outer", suffixes=("_R", "_py"), indicator=True
    )
    n_ref_only = int((m["_merge"] == "left_only").sum())
    n_py_only = int((m["_merge"] == "right_only").sum())
    both = m[m["_merge"] == "both"].copy()

    result = {"n_reaches_R": len(ref), "n_reaches_py": len(py), "n_matched": len(both),
              "reaches_only_in_R": n_ref_only, "reaches_only_in_py": n_py_only}
    worst = {}
    n_over = 0
    for col in COMPARE_COLS:
        a = both[f"{col}_R"].to_numpy(dtype=float)
        b = both[f"{col}_py"].to_numpy(dtype=float)
        abs_err = np.abs(a - b)
        # treat (NaN, NaN) as a match (both models legitimately return NaN for some reaches)
        matched_nan = _both_nan(a, b)
        abs_err = np.where(matched_nan, 0.0, abs_err)
        # relative error against the R reference magnitude (guard tiny denominators)
        denom = np.maximum(np.abs(a), 1e-12)
        rel_err = np.where(matched_nan, 0.0, abs_err / denom)
        # a reach fails if it is over BOTH abs and rel tol, or exactly one side is NaN
        one_side_nan = np.isnan(a) ^ np.isnan(b)
        fail_mask = ((abs_err > abs_tol) & (rel_err > rel_tol)) | one_side_nan
        n_over += int(np.nansum(fail_mask))
        result[f"{col}_max_abs_err"] = float(np.nanmax(abs_err)) if len(abs_err) else np.nan
        result[f"{col}_median_abs_err"] = float(np.nanmedian(abs_err)) if len(abs_err) else np.nan
        result[f"{col}_max_rel_err"] = float(np.nanmax(rel_err)) if len(rel_err) else np.nan
        worst[col] = int(np.nansum(fail_mask))

    result["n_reaches_over_tol"] = n_over
    # require at least one matched reach: two header-only (empty) outputs must NOT pass silently.
    result["passed"] = bool(len(both) > 0 and n_over == 0 and n_ref_only == 0 and n_py_only == 0)
    result["fails_by_col"] = worst
    return result


def main():
    ap = argparse.ArgumentParser(description="Compare R reference vs Python-port model output.")
    ap.add_argument("--ref-dir", required=True, help="dir with <huc4>_R_output.csv")
    ap.add_argument("--py-dir", required=True, help="dir with <huc4>_py_output.csv")
    ap.add_argument("--huc4", nargs="+", required=True, help="basins to compare, e.g. 0107 0108")
    ap.add_argument("--abs-tol", type=float, default=1e-6, help="absolute tolerance [ppm]")
    ap.add_argument("--rel-tol", type=float, default=1e-6, help="relative tolerance")
    ap.add_argument("--report", default=None, help="optional CSV path for the per-basin report")
    args = ap.parse_args()

    ref_dir, py_dir = Path(args.ref_dir), Path(args.py_dir)
    rows = []
    all_pass = True
    print(f"{'HUC4':>8} {'matched':>8} {'CO2 max|err|':>13} {'CO2 max rel':>12} {'over tol':>9} {'result':>7}")
    for huc4 in args.huc4:
        ref_csv = ref_dir / f"{huc4}_R_output.csv"
        py_csv = py_dir / f"{huc4}_py_output.csv"
        if not ref_csv.exists() or not py_csv.exists():
            print(f"{huc4:>8}   MISSING  ({'R' if not ref_csv.exists() else ''}"
                  f"{' py' if not py_csv.exists() else ''})")
            all_pass = False
            rows.append({"huc4": huc4, "passed": False, "note": "missing output file"})
            continue
        r = compare_basin(ref_csv, py_csv, args.abs_tol, args.rel_tol)
        r["huc4"] = huc4
        rows.append(r)
        all_pass = all_pass and r["passed"]
        print(f"{huc4:>8} {r['n_matched']:>8} {r['CO2_ppm_max_abs_err']:>13.3e} "
              f"{r['CO2_ppm_max_rel_err']:>12.3e} {r['n_reaches_over_tol']:>9} "
              f"{'PASS' if r['passed'] else 'FAIL':>7}")

    if args.report:
        pd.DataFrame(rows).to_csv(args.report, index=False)
        print(f"\nper-basin report written to {args.report}")

    print(f"\nOVERALL: {'ALL BASINS MATCH R within tolerance' if all_pass else 'MISMATCH(ES) FOUND'}"
          f"  (abs_tol={args.abs_tol}, rel_tol={args.rel_tol})")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
