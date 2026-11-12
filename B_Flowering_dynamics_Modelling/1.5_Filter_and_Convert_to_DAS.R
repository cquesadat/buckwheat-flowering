#!/usr/bin/env Rscript
# Apply recorded plot-level quality control and convert event dates to days after sowing.
# The fitted models and their parameters are not modified by this script.

input_dir <- "flowering_parameters"
data_dir <- "data"
threshold <- 0.5

read_sowing_dates <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("Trial", "Sowing_Date") %in% names(df))) {
    stop("Sowing-date table must contain Trial and Sowing_Date")
  }
  df$Sowing_Date <- parse_dates(df$Sowing_Date, "Sowing_Date")
  if (anyDuplicated(df$Trial)) stop("Duplicate experiment in sowing-date table")
  df
}

parse_dates <- function(values, field) {
  x <- as.character(values)
  out <- as.Date(rep(NA_character_, length(x)))
  allowed <- c("%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y", "%d-%m-%Y")
  for (fmt in allowed) {
    ix <- which(is.na(out) & !is.na(x) & nzchar(x))
    if (length(ix)) out[ix] <- as.Date(x[ix], format = fmt)
  }
  invalid <- !is.na(x) & nzchar(x) & is.na(out)
  if (any(invalid)) stop("Invalid dates in ", field, ": ", paste(head(unique(x[invalid]), 4), collapse = ", "))
  out
}

mask_abscission <- function(df, bad) {
  # Linear-fit outcomes and durations involving the projected end of abscission.
  cols <- grep("^abscission_|_to_end_abscission$", names(df), value = TRUE)
  if (length(cols)) df[bad, cols] <- NA
  df
}

manual_rules <- read.csv(file.path(data_dir, "manual_exclusions.csv"), stringsAsFactors = FALSE)
expected <- c("experiment", "scope", "field", "value")
if (!all(expected %in% names(manual_rules))) stop("Invalid manual exclusions file")
if (any(!manual_rules$scope %in% c("all", "abscission"))) stop("Invalid exclusion scope")
if (any(!manual_rules$field %in% c("ID", "Plot"))) stop("Invalid exclusion field")

sowing <- read_sowing_dates(file.path(data_dir, "sowing_dates.csv"))
experiments <- c("ETH23", "ETH24", "LfL23", "LfL24", "ETH24L")
files <- file.path(input_dir, paste0("parameters_", experiments, ".csv"))
files <- files[file.exists(files)]
if (!length(files)) stop("No fitted parameter files in ", input_dir)
summary_rows <- list()

for (file in files) {
  exp <- sub("^parameters_(.*)[.]csv$", "\\1", basename(file))
  metrics_path <- file.path(input_dir, paste0("metrics_", exp, ".csv"))
  if (!file.exists(metrics_path)) stop("Missing fit metrics: ", metrics_path)
  params <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
  key <- c("Plot", "ID")
  required_metrics <- c(key, "logistic_R2", "linear_R2")
  if (!all(key %in% names(params)) || !all(required_metrics %in% names(metrics))) {
    stop("Missing Plot/ID keys or R2 metrics for ", exp)
  }
  if (anyDuplicated(params[key]) || anyDuplicated(metrics[key])) {
    stop("Non-unique Plot/ID keys in ", exp)
  }
  k1 <- paste(params$Plot, params$ID, sep = "\r")
  k2 <- paste(metrics$Plot, metrics$ID, sep = "\r")
  ix <- match(k1, k2)
  if (anyNA(ix)) stop("Missing model-fit metrics for some plots in ", exp)
  r2_log <- metrics$logistic_R2[ix]
  r2_lin <- metrics$linear_R2[ix]
  qc <- params[key]
  qc$logistic_R2 <- r2_log
  qc$linear_R2 <- r2_lin
  qc$failed_logistic <- is.na(r2_log) | r2_log < threshold
  qc$failed_abscission <- is.na(r2_lin) | r2_lin < threshold
  qc$manual_remove <- FALSE
  qc$manual_abscission <- FALSE
  qc$late_abscission <- FALSE
  qc$reason <- ""

  local_rules <- manual_rules[manual_rules$experiment == exp, , drop = FALSE]
  for (i in seq_len(nrow(local_rules))) {
    rule <- local_rules[i, ]
    hit <- as.character(params[[rule$field]]) == as.character(rule$value)
    hit[is.na(hit)] <- FALSE
    if (rule$scope == "all") qc$manual_remove <- qc$manual_remove | hit
    if (rule$scope == "abscission") qc$manual_abscission <- qc$manual_abscission | hit
    qc$reason[hit] <- paste0(qc$reason[hit], ifelse(nzchar(qc$reason[hit]), "; ", ""),
                             rule$scope, ":", rule$field, "=", rule$value)
  }
  qc$excluded <- qc$failed_logistic | qc$manual_remove

  # Historical LfL criterion: the projected end of abscission must not pass
  # day of year 335. Apply only to the abscission outcomes, after fitting.
  if (startsWith(exp, "LfL") && "abscission_projected_end_date" %in% names(params)) {
    projected <- parse_dates(params$abscission_projected_end_date, "abscission_projected_end_date")
    cutoff <- as.integer(format(projected, "%j"))
    qc$late_abscission <- !is.na(cutoff) & cutoff > 335L
  }
  bad_abscission <- qc$failed_abscission | qc$manual_abscission | qc$late_abscission
  params <- mask_abscission(params, bad_abscission)
  params <- params[!qc$excluded, , drop = FALSE]

  # Event timing variables are expressed as DAS in the GWAS input tables.
  sowing_date <- sowing$Sowing_Date[match(exp, sowing$Trial)]
  if (length(sowing_date) != 1L || is.na(sowing_date)) stop("Missing sowing date for ", exp)
  date_cols <- grep("_date$", names(params), value = TRUE)
  for (date_col in date_cols) {
    event <- parse_dates(params[[date_col]], date_col)
    params[[sub("_date$", "_das", date_col)]] <- as.integer(event - sowing_date)
  }
  params[date_cols] <- NULL
  clean_path <- file.path(input_dir, paste0("parameters_", exp, "_clean_das.csv"))
  # Keep only identifiers and extracted numeric traits for downstream SpATS.
  params <- params[, c("Plot", "ID", setdiff(names(params), c("Plot", "ID"))), drop = FALSE]
  write.csv(params, clean_path, row.names = FALSE, na = "")
  write.csv(qc, file.path(input_dir, paste0("quality_control_", exp, ".csv")), row.names = FALSE)
  summary_rows[[exp]] <- data.frame(experiment = exp, n_fit = nrow(qc),
    n_low_logistic_R2 = sum(qc$failed_logistic),
    n_low_linear_R2 = sum(qc$failed_abscission),
    n_manual_excluded = sum(qc$manual_remove),
    n_manual_abscission = sum(qc$manual_abscission),
    n_late_abscission = sum(qc$late_abscission),
    n_remaining = nrow(params))
  message(exp, ": retained ", nrow(params), " of ", nrow(qc), " fitted plots")
}
write.csv(do.call(rbind, summary_rows), file.path(input_dir, "quality_control_summary.csv"), row.names = FALSE)
