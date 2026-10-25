#!/usr/bin/env python
"""
Phase B: compare BOTH the R reference output and the Python port output against
Craig's saved per-basin result (cache/final_<huc4>.csv), reach by reach.

Phase A (compare.py) proves R == Python on identical inputs. This asks the stronger
question: does each reproduce Craig's ACTUAL published output? final_<huc4>.csv was
produced with the CALIBRATED params, so run export_reference.R with
CONUS_PARAM_MODE=calib (R and Python then use the same params Craig used).

Interpretation
  - R  vs final small  -> our invocation (data, AHG, upstreamDF=NA) reproduces Craig's pipeline.
  - py vs final small  -> the port reproduces the real result.
  - If R vs final diverges while R == py (Phase A), the gap is an INPUT/import difference
    (e.g. cross-basin CO2 import, which a single-basin run sets to NA), NOT a port bug.

"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

FIELDS = ["CO2_ppm", "FCO2_gC_m2_yr"]


def _load(path, cols):
    df = pd.read_csv(path)
    if "NHDPlusID" not in df.columns:
        raise SystemExit(f"{path}: no NHDPlusID column")
    df["NHDPlusID"] = pd.to_numeric(df["NHDPlusID"]).astype("int64")
    keep = ["NHDPlusID"] + [c for c in cols if c in df.columns]
    return df[keep]


def _err_stats(pred, ref, abs_tol, rel_tol):
    pred = np.asarray(pred, dtype=float)
    ref = np.asarray(ref, dtype=float)
    ok = np.isfinite(pred) & np.isfinite(ref)
    absd = np.abs(pred - ref)[ok]
    reld = (np.abs(pred - ref) / np.maximum(np.abs(ref), 1e-12))[ok]
    n_over = int(((absd > abs_tol) & (reld > rel_tol)).sum())
    return (
        float(absd.max()) if absd.size else float("nan"),
        float(reld.max()) if reld.size else float("nan"),
        n_over,
        int((~ok).sum()),
    )


def compare_basin(huc4, ref_dir, final_dir, abs_tol, rel_tol):
    r = _load(os.path.join(ref_dir, f"{huc4}_R_output.csv"), FIELDS)
    p = _load(os.path.join(ref_dir, f"{huc4}_py_output.csv"), FIELDS)
    f = _load(os.path.join(final_dir, f"final_{huc4}.csv"), FIELDS)

    row = {"huc4": huc4, "n_R": len(r), "n_py": len(p), "n_final": len(f)}
    m = r.merge(p, on="NHDPlusID", suffixes=("_R", "_py")).merge(f, on="NHDPlusID")
    row["n_matched"] = len(m)
    for fld in FIELDS:
        for tag in ("R", "py"):
            a, rl, over, nbad = _err_stats(m[f"{fld}_{tag}"], m[fld], abs_tol, rel_tol)
            row[f"{tag}_vs_final.{fld}.max_abs"] = a
            row[f"{tag}_vs_final.{fld}.max_rel"] = rl
            row[f"{tag}_vs_final.{fld}.n_over"] = over
        # Phase A recap (R vs py) on the same matched set
        a, rl, over, nbad = _err_stats(m[f"{fld}_R"], m[f"{fld}_py"], abs_tol, rel_tol)
        row[f"R_vs_py.{fld}.max_abs"] = a
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description="Phase B: R and Python vs Craig's final_<huc4>.csv")
    ap.add_argument("--ref-dir", required=True, help="dir with <huc4>_R_output.csv and <huc4>_py_output.csv")
    ap.add_argument("--final-dir", required=True, help="dir with Craig's final_<huc4>.csv")
    ap.add_argument("--huc4", nargs="+", required=True)
    ap.add_argument("--abs-tol", type=float, default=1e-6)
    ap.add_argument("--rel-tol", type=float, default=1e-6)
    ap.add_argument("--report", default=None, help="optional CSV report path")
    args = ap.parse_args(argv)

    rows = []
    for huc4 in args.huc4:
        rows.append(compare_basin(huc4, args.ref_dir, args.final_dir, args.abs_tol, args.rel_tol))

    rep = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 60)

    for r in rows:
        print(f"\n===== basin {r['huc4']} =====")
        print(f"  reaches: R={r['n_R']}  py={r['n_py']}  final={r['n_final']}  matched(3-way)={r['n_matched']}")
        for fld in FIELDS:
            print(f"  {fld}:")
            print(f"    R  vs final : max|abs|={r[f'R_vs_final.{fld}.max_abs']:.6g}"
                  f"  max|rel|={r[f'R_vs_final.{fld}.max_rel']:.6g}"
                  f"  n_over_tol={r[f'R_vs_final.{fld}.n_over']}")
            print(f"    py vs final : max|abs|={r[f'py_vs_final.{fld}.max_abs']:.6g}"
                  f"  max|rel|={r[f'py_vs_final.{fld}.max_rel']:.6g}"
                  f"  n_over_tol={r[f'py_vs_final.{fld}.n_over']}")
            print(f"    R  vs py    : max|abs|={r[f'R_vs_py.{fld}.max_abs']:.6g}   (Phase A recap)")

    if args.report:
        os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
        rep.to_csv(args.report, index=False)
        print(f"\nwrote {args.report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
