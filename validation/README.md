# validation

These scripts check that the Python port gives the same reach-by-reach CO2 and emissions as the original R model when both get the same prepared hydrography and parameters.

There are three steps, submitted as a Slurm dependency chain by `run_all.sh` (run it from the repo root):

1. `slurm/1_export_reference.sbatch` runs the original R model through `R/export_reference.R` and writes reference CSVs for each basin.
2. `slurm/2_run_python.sbatch` runs the port on the same inputs with `python/run_python_model.py`.
3. `slurm/3_compare.sbatch` diffs the two with `python/compare.py` and exits non-zero if any basin does not match.

The report ends up in `out/compare_report.csv` and the logs in `logs/`. `slurm/run_validation_container.sbatch` does all three in one job, with R in a rocker container.

Basins, paths and tolerances are set in `config.sh`. I ran 0107, 0308, 0505, 0704, 0805, 1104, 1703 and 1806.
