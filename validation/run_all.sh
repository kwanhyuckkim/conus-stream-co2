#!/usr/bin/env bash
# Submit the 3-step validation as a SLURM dependency chain: export(R) -> run(Python) -> compare.
# Run from the repo root on Unity:  bash validation/run_all.sh
set -euo pipefail
source validation/config.sh
mkdir -p validation/logs validation/out

j1=$(sbatch --parsable validation/slurm/1_export_reference.sbatch)
echo "submitted export(R)     : $j1"
j2=$(sbatch --parsable --dependency=afterok:$j1 validation/slurm/2_run_python.sbatch)
echo "submitted run(Python)   : $j2  (after $j1)"
j3=$(sbatch --parsable --dependency=afterok:$j2 validation/slurm/3_compare.sbatch)
echo "submitted compare       : $j3  (after $j2)"
echo
echo "watch:   squeue -u \$USER"
echo "result:  validation/logs/compare_${j3}.out  and  validation/out/compare_report.csv"
