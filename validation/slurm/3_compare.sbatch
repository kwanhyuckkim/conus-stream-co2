#!/usr/bin/env bash
#SBATCH --job-name=conus_compare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=validation/logs/compare_%j.out
#SBATCH --error=validation/logs/compare_%j.err
#
# Step 3/3: diff R vs Python reach by reach. Exits non-zero if any basin mismatches.

set -euo pipefail
source validation/config.sh

module load conda/latest 2>/dev/null || module load miniconda 2>/dev/null || true
eval "$(conda shell.bash hook)"   # init conda for non-interactive shells (adjust module name for Unity)
conda activate "$PY_ENV"

python validation/python/compare.py \
  --ref-dir "$VALIDATION_OUT_R" \
  --py-dir  "$VALIDATION_OUT_PY" \
  --huc4    $CONUS_BASINS \
  --abs-tol "$ABS_TOL" --rel-tol "$REL_TOL" \
  --report  "validation/out/compare_report.csv"
