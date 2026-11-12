#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: 0.1_Add_DAS.R <path-to-Results_Flowering_Perc_<EXP>.csv>")
}
infile <- args[1]
if (!file.exists(infile)) stop("Input file not found: ", infile)

basedir <- dirname(infile)
sowing_file <- file.path(basedir, "sowing_dates.csv")
if (!file.exists(sowing_file)) stop("Sowing dates file not found in data/: ", sowing_file)

df <- read.csv(infile, stringsAsFactors = FALSE)
sow <- read.csv(sowing_file, stringsAsFactors = FALSE)

parse_dates_vec <- function(v) {
  res <- as.Date(rep(NA, length(v)))
  formats <- c("%d.%m.%Y", "%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y")
  for (fmt in formats) {
    nas <- is.na(res)
    if (!any(nas)) break
    parsed <- as.Date(v[nas], format = fmt)
    res[nas] <- parsed
  }
  if (any(is.na(res)) && requireNamespace("lubridate", quietly = TRUE)) {
    nas <- is.na(res)
    res[nas] <- as.Date(lubridate::dmy(v[nas]))
  }
  return(res)
}

if (!"Date" %in% colnames(df)) stop("Input CSV must have a 'Date' column")
df$Date_parsed <- parse_dates_vec(df$Date)

# derive experiment name from filename: Results_Flowering_Perc_<EXP>.csv
expname <- sub('^Results_Flowering_Perc_(.*)\\.csv$', '\\1', basename(infile))
srow <- sow[sow$Trial == expname, , drop = FALSE]

if (nrow(srow) == 0) {
  warning("No sowing entry found for experiment: ", expname, " — DAS set to NA")
  df$DAS <- NA
} else {
  sow_date_str <- srow$Sowing_Date[1]
  sow_date <- as.Date(NA)
  # try common formats
  try_formats <- c("%d.%m.%Y", "%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y")
  for (fmt in try_formats) {
    if (is.na(sow_date)) sow_date <- as.Date(sow_date_str, format = fmt)
  }
  if (is.na(sow_date) && requireNamespace("lubridate", quietly = TRUE)) {
    sow_date <- as.Date(lubridate::dmy(sow_date_str))
  }
  if (is.na(sow_date)) {
    warning("Could not parse sowing date for experiment: ", expname, " — DAS set to NA")
    df$DAS <- NA
  } else {
    # Use difference in days between image date and sowing date
    df$DAS <- as.integer(difftime(df$Date_parsed, sow_date, units = "days"))
  }
}

# remove helper column and write back
df$Date_parsed <- NULL
write.csv(df, infile, row.names = FALSE)
cat("Wrote DAS into:", infile, "\n")
