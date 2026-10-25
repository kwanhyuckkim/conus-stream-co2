# conus-stream-co2

A Python port of Craig Brinkerhoff's CONUS_carbon R model, which simulates CO2 in river networks across the contiguous US.

- Original R code: https://github.com/craigbrinkerhoff/CONUS_carbon
- Method paper: Saccardi, Brinkerhoff, et al. (2024), Modeling CO2 dynamics in river networks across the contiguous United States, AGU Advances, https://doi.org/10.1029/2024AV001294

All scientific credit for the model goes to those authors. I ported the forward model to Python so I could run parameter sweeps and calibration experiments for my postdoc work with Matt Winnick at UMass Amherst. The calibration runs and their results are not in this repo.

## What's here

- `conus_carbon/`: the model (steady-state pCO2 for each reach, solved from upstream to downstream over a basin), the calibration cost, hydraulics, I/O, a parameter sweep, and a small synthetic basin for testing
- `tests/`: physics sanity checks and a regression test
- `validation/`: the scripts I used to run the R original and this port on the same hydrography and parameters and compare them reach by reach
- `CONUS_carbon_flow.ipynb`: the whole model run on the synthetic basin, no lab data needed
- `CONUS_carbon_sweep.ipynb`: a linearity and influence sweep on that basin

## How close it is to the R code

On the basins I checked, the Python and R outputs agree to about 1e-7 ppm. The code changes are only for speed, mainly an upstream-index map that is built once instead of filtering the reach table for every reach. Two R behaviors are copied on purpose so the two codes match: the lake k600 function gets lake area in m2 while its bins are in km2, and the biozone term (k_bz) adds CO2 as a one-way source. The model authors confirmed the second one is intended.

`setup_hydrography` (reading the NHD-HR geodatabases) is not ported. Both models are fed the hydrography table that the R code prepares.

## Running it

```bash
pip install -e .
pytest tests/
```

The notebooks need `pip install -e ".[notebook]"`. The validation scripts expect the R original and the NHD-HR data on our cluster, so check the paths in `validation/config.sh` first.
