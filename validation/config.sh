#!/usr/bin/env bash
# Shared config for the R-vs-Python validation harness. Sourced by run_all.sh and the .sbatch files.

# Space-separated HUC4 test basins. An EXTERNAL value wins, so you can run one basin per job:
#   sbatch --export=ALL,CONUS_BASINS="0308" validation/slurm/run_validation_container.sbatch
# Validated 2026-07: 0107 0308 0505 0704 0805 1104 1703 1806 (Indiana basins + 1710 need extra files).
export CONUS_BASINS="${CONUS_BASINS:-0107}"

# NHD-HR data root on Unity. Layout must be {ROOT}/HUC2_{huc2}/NHDPLUS_H_{huc4}_HU4_GDB/....gdb
# CONFIRMED on Unity 2026-07-14: the shared lab copy has HUC2_01..19 and basin 0107's gdb.
export CONUS_DATA_ROOT="/nas/cee-ice/data/NHDPlus_HR"

# Path to the ORIGINAL R model source (Craig's src/ with model.R/utils.R/...). This repo ships the
# Python port + validation only, NOT Craig's R model. CONFIRMED on Unity 2026-07-14: Craig's src is here.
export CONUS_SRC="/nas/cee-water/cjgleason/craig/CONUS_carbon/src"

# AHG dead-path fix. src/setup_hydrography.R:224-225 hardcode depAHG.rds/widAHG.rds under
#   /nas/cee-water/cjgleason/craig/RSK600/cache/  -> THAT PATH NO LONGER EXISTS, so setupHydrography
#   crashes. The same RSK600 .rds live under heejin's copy (CONUS_AHG_DIR below). To run the R side,
#   copy Craig's src to a WRITABLE dir and repoint those two lines:
#     cp -r "$CONUS_SRC" ./R_src && export CONUS_SRC="$PWD/R_src"
#     sed -i "s#/nas/cee-water/cjgleason/craig/RSK600/cache#$CONUS_AHG_DIR#g" R_src/setup_hydrography.R
export CONUS_AHG_DIR="/nas/cee-water/cjgleason/heejin/CONUS_carbon_CT/RSK600/cache"

# Where the standardized CSVs go
export VALIDATION_OUT_R="validation/out/R"
export VALIDATION_OUT_PY="validation/out/python"

# Comparison tolerances (a faithful port should match R to ~1e-9; 1e-6 leaves margin for CSV round-trip).
export ABS_TOL="1e-6"
export REL_TOL="1e-6"

# Conda env names on Unity (used only by the generic conda .sbatch in slurm/1_,2_,3_).
export R_ENV="env-CONUS_carbon"        # the R pipeline env (environment.yml)
export PY_ENV="conus_carbon_py"        # a python env with numpy + pandas (see requirements.txt)

# --- What actually worked on Unity 2026-07-16 (used by slurm/run_validation_container.sbatch) ---
# The bare r/4.4.0 module has NO sf; the rocker ml-verse apptainer container does (sf 1.0.16,
# dplyr, readr). Bind /nas so the containerized R can read the NHD-HR data. Python runs from a
# direct env path (no conda activate needed). Adjust CONUS_PY_BIN to any env with numpy+pandas.
export CONUS_R_SIF="/modules/apps/ood/rstudio-server/ml-verse_4.4.0.sif"
export CONUS_APPTAINER_BIND="/nas,/work,/home"
export CONUS_PY_BIN="/work/kwanhyuckkim_umass_edu/.conda/envs/swot-urban/bin/python"
export CONUS_PARAM_MODE="calib"        # "calib" -> Craig's params (enables Phase B); "fixed" -> test vector
export CONUS_CALIB_CSV="/nas/cee-water/cjgleason/craig/CONUS_carbon/cache/calibParams.csv"

# NOTE: the SLURM partition is set as `#SBATCH --partition=cpu` in validation/slurm/*.sbatch (a
# #SBATCH directive can't read this file). Confirm it with `sinfo -s` on Unity and edit the three
# .sbatch files if your partition is named differently (e.g. cpu-preempt).
# CONUS_DATA_ROOT / CONUS_SRC / CONUS_AHG_DIR are lab-internal (Unity/NAS): adjust off-cluster.
#
# REFERENCE OUTPUTS found on Unity 2026-07-14 (an independent way to validate the port):
#   /nas/cee-water/cjgleason/craig/CONUS_carbon/cache/final_<huc4>.csv   (211 basins incl. 0107)
#     columns: NHDPlusID, ..., CO2_ppm, ..., FCO2_gC_m2_yr  = Craig's published per-reach outputs.
#   /nas/cee-water/cjgleason/craig/CONUS_carbon/cache/calibParams.csv    = calibrated params per basin
#     (0107: Cbz_riv=267.56, Cbz_lake=0, Fwc_riv=0, Fwc_lake=0).
#   These check Python (run with the calibrated params) against Craig's ACTUAL results, not just
#   R-vs-Python self-consistency. final_*.csv is an OUTPUT subset (no HRT/k_bz/HydroSeq), so it is a
#   comparison TARGET only; the full 24-col hydrography still comes from R (_targets or setupHydrography).
