#############################################################################
## export_reference.R
## Run the ORIGINAL R model on test basins and export standardized CSVs that the
## Python port is validated against, reach by reach.
##
## Run this on Unity (where the NHD-HR data + AHG models + R env live).
##
## Per basin it writes THREE files into --out (default validation/out/R):
##   <huc4>_hydrography.csv : the 24-column prepared hydrography (INPUT to both models)
##   <huc4>_params.csv      : the parameter vector used (Cbz_riv,Cbz_lake,Fwc_riv,Fwc_lake)
##   <huc4>_R_output.csv    : R runModel per-reach output (NHDPlusID, CO2_ppm, FCO2_gC_m2_yr)
##
## The Python side then runs conus_carbon on the SAME hydrography + params and
## validation/python/compare.py diffs the two outputs.
#############################################################################

## ---- CONFIG (edit these)
## basins come from the CONUS_BASINS env var (space-separated) if set, else this default.
basins        <- if (nzchar(Sys.getenv("CONUS_BASINS")))
                   strsplit(trimws(Sys.getenv("CONUS_BASINS")), "\\s+")[[1]] else c("0107")
path_to_data  <- Sys.getenv("CONUS_DATA_ROOT",
                            "/nas/cee-water/cjgleason/craig/CONUS_ephemeral_data")  # NHD-HR root
out_dir       <- Sys.getenv("VALIDATION_OUT", "validation/out/R")
src_dir       <- Sys.getenv("CONUS_SRC", "src")   # folder with model.R/utils.R/... relative to repo root

## Physical constants (must match constants.py)
C_groundwater <- 16000    # [ppm]  (_targets.R:40)
C_atmosphere  <- 400      # [ppm]  (_targets.R:38)
emergenceQ    <- 0.000105 # [m3/s] (_targets.R:39)

## Parameter source (env-driven)
##   CONUS_PARAM_MODE = "fixed" -> use the test vector below (validates forward-model LOGIC only)
##   CONUS_PARAM_MODE = "calib" -> read Craig's calibrated params from CONUS_CALIB_CSV (columns
##       basin,Cbz_riv,Cbz_lake,Fwc_riv,Fwc_lake), so the output can be checked against his saved
##       final_<huc4>.csv reference.
param_mode    <- Sys.getenv("CONUS_PARAM_MODE", "fixed")
calib_csv     <- Sys.getenv("CONUS_CALIB_CSV", "")
fixed_params  <- list(Cbz_riv = 1000.0, Cbz_lake = 1000.0, Fwc_riv = 1e-4, Fwc_lake = 1e-5)
## ---------------------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(src_dir, "utils.R"))
source(file.path(src_dir, "model.R"))
source(file.path(src_dir, "calibrate_model.R"))
source(file.path(src_dir, "run_model.R"))
source(file.path(src_dir, "setup_hydrography.R"))

## Attach the packages Craig's functions expect AT CALL TIME. setup_hydrography.R uses dplyr's
## `%>%` pipe un-namespaced (and sf / readr), so these must be attached, not merely installed
## (matches the @import sf / dplyr / readr roxygen tags). The rocker container provides all three.
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
})

## write a data frame at high precision so the CSV round-trip does not perturb the comparison
write_hi <- function(df, path) {
  old <- options(digits = 17, scipen = 999); on.exit(options(old))
  utils::write.csv(df, path, row.names = FALSE)
}

for (huc4 in basins) {
  message("=== ", huc4, " ===")

  ## 1) prepared hydrography (the 24-column routing table)
  ##    Default: build it from NHD-HR exactly as the pipeline does.
  ##    Alternatives if the data is not mounted: tar_read(<hydrography target>) or read a saved RDS.
  ## 1710a/1710b are SPLITS of basin 1710; the .gdb on disk is named 1710, so build the hydrography
  ## with the BASE id, then split. (Calling setupHydrography("1710a") would look for a non-existent
  ## NHDPLUS_H_1710a_HU4_GDB and crash.)
  base_huc4 <- if (huc4 %in% c("1710a", "1710b")) "1710" else huc4
  hydrography <- setupHydrography(path_to_data, base_huc4)
  if (huc4 %in% c("1710a", "1710b")) hydrography <- split1710(hydrography, huc4)

  ## 2) parameters
  if (param_mode == "calib") {
    if (!nzchar(calib_csv)) stop("CONUS_PARAM_MODE=calib requires CONUS_CALIB_CSV (path to calibParams.csv)")
    cp  <- utils::read.csv(calib_csv, colClasses = "character")   # keep '0107' as a string, not integer 107
    row <- cp[cp$basin == huc4, ]
    if (nrow(row) != 1) stop("calibParams.csv: expected one row for basin ", huc4, ", got ", nrow(row))
    params <- list(Cbz_riv = as.numeric(row$Cbz_riv), Cbz_lake = as.numeric(row$Cbz_lake),
                   Fwc_riv = as.numeric(row$Fwc_riv), Fwc_lake = as.numeric(row$Fwc_lake))
    message("  calib params: Cbz_riv=", params$Cbz_riv, " Cbz_lake=", params$Cbz_lake,
            " Fwc_riv=", params$Fwc_riv, " Fwc_lake=", params$Fwc_lake)
  } else if (param_mode == "fixed") {
    params <- fixed_params
  } else {
    stop("unknown CONUS_PARAM_MODE: ", param_mode, " (use 'fixed' or 'calib')")
  }
  calibrationResults <- list(Cbz_riv = params$Cbz_riv, Cbz_lake = params$Cbz_lake,
                             Fwc_riv = params$Fwc_riv, Fwc_lake = params$Fwc_lake)

  ## 3) run the ORIGINAL R model (forward CO2 + emissions)
  R_out <- runModel(hydrography, calibrationResults, C_groundwater, C_atmosphere,
                    huc4, emergenceQ, NA)

  ## 4) write the three standardized files
  write_hi(as.data.frame(hydrography), file.path(out_dir, paste0(huc4, "_hydrography.csv")))
  write_hi(data.frame(Cbz_riv = params$Cbz_riv, Cbz_lake = params$Cbz_lake,
                      Fwc_riv = params$Fwc_riv, Fwc_lake = params$Fwc_lake,
                      Cgw = C_groundwater, Catm = C_atmosphere, emergenceQ = emergenceQ),
           file.path(out_dir, paste0(huc4, "_params.csv")))
  ref <- as.data.frame(R_out)[, c("NHDPlusID", "CO2_ppm", "FCO2_gC_m2_yr")]
  write_hi(ref, file.path(out_dir, paste0(huc4, "_R_output.csv")))

  message("  wrote hydrography (", nrow(hydrography), " reaches), params, and R output.")
}
message("Done. Reference CSVs in: ", out_dir)
