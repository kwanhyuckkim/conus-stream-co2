#!/bin/bash
#SBATCH --job-name=conus_valid
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=validation/logs/valid_%j.log
#SBATCH --error=validation/logs/valid_%j.log
#
# One-shot R<->Python validation on Unity using the rocker apptainer container for R and a
# direct python interpreter for the port. Verified on basin 0107 (2026-07-16): Phase A PASS,
# Phase B (vs Craig's final_<huc4>.csv) matched to ~1e-8 ppm.
#
# Prereqs
#   - run from the repo ROOT:  sbatch validation/slurm/run_validation_container.sbatch
#   - a writable copy of Craig's R src with the AHG path patched (see validation/config.sh
#     CONUS_AHG_DIR note); set CONUS_SRC to it.
#   - the CONUS boundary shapefile copied into ./data/ (setup_hydrography.R reads data/cb_2018_*).
#
set -uo pipefail
source validation/config.sh

OUTR="$VALIDATION_OUT_R"
OUT="$(dirname "$OUTR")"
mkdir -p "$OUTR" validation/logs
export VALIDATION_OUT="$OUTR"
FINAL_DIR="$(dirname "$CONUS_CALIB_CSV")"      # Craig's final_<huc4>.csv live next to calibParams.csv
TAG="$(echo "$CONUS_BASINS" | tr ' ' '_')"    # tag reports by basin(s) so concurrent per-basin jobs don't clobber

# --- self-setup (idempotent) -------------------------------------------------------------------
# (a) AHG-patched, WRITABLE copy of the R src. Craig's src is read-only and setup_hydrography.R:224-225
#     hardcode a dead RSK600 path; repoint those two lines to CONUS_AHG_DIR (a live copy of the same .rds).
RSRC="$HOME/conus_validation/R_src"
if [ ! -f "$RSRC/setup_hydrography.R" ] || ! grep -q "$CONUS_AHG_DIR" "$RSRC/setup_hydrography.R"; then
  echo "building AHG-patched R src at $RSRC (from $CONUS_SRC)"
  rm -rf "$RSRC"; cp -r "$CONUS_SRC" "$RSRC"
  sed -i "s#/nas/cee-water/cjgleason/craig/RSK600/cache#$CONUS_AHG_DIR#g" "$RSRC/setup_hydrography.R"
fi
export CONUS_SRC="$RSRC"                       # <- R sources the patched copy, not Craig's read-only src
# (b) CONUS boundary shapefile: setup_hydrography.R reads it by RELATIVE path data/cb_2018_us_nation_20m.shp
[ -f data/cb_2018_us_nation_20m.shp ] || cp "$CONUS_DATA_ROOT"/other_shapefiles/cb_2018_us_nation_20m.* data/ 2>/dev/null || true
# -----------------------------------------------------------------------------------------------

echo "########## STEP 1: R export (rocker apptainer; ${CONUS_PARAM_MODE} params) ##########"; date
apptainer exec --bind "$CONUS_APPTAINER_BIND" "$CONUS_R_SIF" \
  Rscript validation/R/export_reference.R
r1=$?; echo "R export exit=$r1"; ls -l "$OUTR"
[ $r1 -eq 0 ] || { echo "R export failed; stopping."; exit $r1; }

echo "########## STEP 2: Python model on the SAME hydrography + params ##########"
for h in $CONUS_BASINS; do
  "$CONUS_PY_BIN" validation/python/run_python_model.py \
    --hydrography "$OUTR/${h}_hydrography.csv" --params "$OUTR/${h}_params.csv" \
    --huc4 "$h" --out "$OUTR/${h}_py_output.csv"
done

echo "########## STEP 3: Phase A  (R vs Python, expect ~1e-9 match) ##########"
"$CONUS_PY_BIN" validation/python/compare.py \
  --ref-dir "$OUTR" --py-dir "$OUTR" --huc4 $CONUS_BASINS \
  --abs-tol "$ABS_TOL" --rel-tol "$REL_TOL" --report "$OUT/compareA_${TAG}.csv" || echo "[Phase A exit nonzero]"

echo "########## STEP 4: Phase B  (R & Python vs Craig final_<huc4>.csv) ##########"
"$CONUS_PY_BIN" validation/python/compare_to_reference.py \
  --ref-dir "$OUTR" --final-dir "$FINAL_DIR" --huc4 $CONUS_BASINS \
  --abs-tol "$ABS_TOL" --rel-tol "$REL_TOL" --report "$OUT/compareB_${TAG}.csv"

echo "########## DONE ##########"; date
