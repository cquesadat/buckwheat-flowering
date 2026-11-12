#!/usr/bin/env Rscript
# Fit flowering and abscission models for each plot in each experiment.

library(lubridate)
library(dplyr)

data_dir <- file.path("data")
out_dir <- file.path("flowering_parameters")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("utils", "FUN_model_fitting.R"))

files <- list.files(data_dir, pattern = "^Results_Flowering_Perc_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  message("No Results_Flowering_Perc_*.csv files found in data/")
  quit(status = 0)
}

for (f in files) {
  message("Processing: ", f)
  df <- read.csv(f, stringsAsFactors = FALSE)
  if (!all(c("Plot", "ID", "Date") %in% names(df))) {
    stop("Required columns Plot, ID and Date are missing from: ", f)
  }

  if (!"Total_Area" %in% names(df)) {
    stop("Required column Total_Area is missing from: ", f)
  }

  raw_dates <- as.character(df$Date)
  parsed_dates <- as.Date(rep(NA_character_, length(raw_dates)))
  for (fmt in c("%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y", "%d-%m-%Y")) {
    missing <- is.na(parsed_dates)
    parsed_dates[missing] <- as.Date(raw_dates[missing], format = fmt)
  }
  if (any(is.na(parsed_dates))) stop("Unrecognized or missing dates in: ", f)
  df$Date <- parsed_dates

  # experiment name from filename
  expname <- sub('^Results_Flowering_Perc_(.*)\\.csv$', '\\1', basename(f))

  # Manual exclusions and fit-quality thresholds are applied to fitted traits in
  # 1.5_Filter_and_Convert_to_DAS.R, not to the original measurements.

  # build genotype list keyed by ID
  genotype_list <- split(df, df$ID)

  # date range for fitting
  date_range <- range(df$Date, na.rm = TRUE)

  # The model-fitting function sources its helper functions from utils/.
  res_list <- lapply(genotype_list, function(g) {
    tryCatch({
      get_parameters_start_to_max(g, date_range)
    }, error = function(e) {
      message("Error fitting genotype ", unique(g$ID), ": ", e$message)
      return(NULL)
    })
  })

  # collect parameters and metrics
  params <- do.call("rbind", lapply(res_list, function(x) if(!is.null(x)) x$parameters else NULL))
  metrics <- do.call("rbind", lapply(res_list, function(x) if(!is.null(x) && !is.null(x$metrics)) x$metrics else NULL))

  out_params <- file.path(out_dir, paste0("parameters_", expname, ".csv"))
  out_metrics <- file.path(out_dir, paste0("metrics_", expname, ".csv"))
  if (!is.null(params) && nrow(params) > 0) write.csv(params, out_params, row.names = FALSE)
  if (!is.null(metrics) && nrow(metrics) > 0) write.csv(metrics, out_metrics, row.names = FALSE)
  message("Wrote:", ifelse(!is.null(params) && nrow(params)>0, out_params, " (no params)"), ", ", ifelse(!is.null(metrics) && nrow(metrics)>0, out_metrics, " (no metrics)"))
}

message("All done. Per-experiment outputs in: ", out_dir)
