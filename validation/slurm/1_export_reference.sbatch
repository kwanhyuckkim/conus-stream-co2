#!/usr/bin/env bash
#SBATCH --job-name=conus_export_R
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=validation/logs/export_R_%j.out
#SBATCH --error=validation/logs/export_R_%j.err
#
# Step 1/3: run the ORIGINAL R model and export reference CSVs. Run from the repo root on Unity.

set -euo pipefail
source validation/config.sh
mkdir -p "$VALIDATION_OUT_R" validation/logs

module load conda/latest 2>/dev/null || module load miniconda 2>/dev/null || true
eval "$(conda shell.bash hook)"   # init conda for non-interactive shells (adjust module name for Unity)
conda activate "$R_ENV"

echo "[export] basins: $CONUS_BASINS"
Rscript validation/R/export_reference.R
echo "[export] done -> $VALIDATION_OUT_R"
