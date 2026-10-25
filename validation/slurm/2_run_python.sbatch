#!/usr/bin/env bash
#SBATCH --job-name=conus_run_py
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=validation/logs/run_py_%j.out
#SBATCH --error=validation/logs/run_py_%j.err
#
# Step 2/3: run the Python port on the SAME hydrography + params R exported. Run from repo root.

set -euo pipefail
source validation/config.sh
mkdir -p "$VALIDATION_OUT_PY" validation/logs

module load conda/latest 2>/dev/null || module load miniconda 2>/dev/null || true
eval "$(conda shell.bash hook)"   # init conda for non-interactive shells (adjust module name for Unity)
conda activate "$PY_ENV"

for huc4 in $CONUS_BASINS; do
  python validation/python/run_python_model.py \
    --hydrography "$VALIDATION_OUT_R/${huc4}_hydrography.csv" \
    --params      "$VALIDATION_OUT_R/${huc4}_params.csv" \
    --huc4        "$huc4" \
    --out         "$VALIDATION_OUT_PY/${huc4}_py_output.csv"
done
echo "[python] done -> $VALIDATION_OUT_PY"
