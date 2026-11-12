#!/usr/bin/env Rscript

# Estimate chromosome-specific, background-adjusted LD-decay windows.
# The chromosome-level LD calculations are in 5.1_Estimate_LD_Decay_Chromosome.R.

suppressPackageStartupMessages(library(data.table))

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_path) != 1L) {
  stop("Run with Rscript 5_Estimate_LD_Decay.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path), winslash = "/"))
project_root <- script_dir

snp_file <- file.path(project_root, "data", "SNP_information.txt")
genotype_file <- file.path(project_root, "data", "GWAFFs_imputed.txt")
worker <- file.path(script_dir, "5.1_Estimate_LD_Decay_Chromosome.R")
out_dir <- file.path(project_root, "outputs", "LD_decay")

for (path in c(snp_file, genotype_file, worker)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

cli <- commandArgs(trailingOnly = TRUE)
if (length(cli) > 1L || (length(cli) == 1L && cli != "--resume")) {
  stop("Usage: Rscript 5_Estimate_LD_Decay.R [--resume]")
}
resume <- length(cli) == 1L

snp_info <- fread(snp_file, data.table = TRUE)
chr_col <- grep("Cromosome|Chromosome|Chr", names(snp_info), ignore.case = TRUE, value = TRUE)[1]
if (is.na(chr_col)) stop("SNP_information.txt is missing a chromosome column")
chromosomes <- sort(unique(suppressWarnings(as.integer(snp_info[[chr_col]]))))
chromosomes <- chromosomes[is.finite(chromosomes)]
if (length(chromosomes) == 0L) stop("No numeric chromosomes found in SNP_information.txt")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) {
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
}
if (!file.exists(rscript)) stop("Cannot locate Rscript")

summaries <- vector("list", length(chromosomes))
binned_files <- character(length(chromosomes))
smooth_files <- character(length(chromosomes))

for (i in seq_along(chromosomes)) {
  chr <- chromosomes[i]
  chr_dir <- file.path(out_dir, paste0("chr", chr))
  summary_file <- file.path(chr_dir, paste0("ld_decay_window_chr", chr, ".csv"))
  binned_files[i] <- file.path(chr_dir, paste0("ld_decay_binned_chr", chr, ".csv"))
  smooth_files[i] <- file.path(chr_dir, paste0("ld_decay_smooth_chr", chr, ".csv"))

  if (!resume || !file.exists(summary_file) || !file.exists(binned_files[i])) {
    message("Estimating LD decay on chromosome ", chr, " (", i, "/", length(chromosomes), ")")
    status <- system2(rscript, c(shQuote(worker), "--chr", as.character(chr)))
    if (status != 0L) stop("LD-decay estimation failed on chromosome ", chr)
  } else {
    message("Reusing LD-decay results for chromosome ", chr)
  }

  if (!file.exists(summary_file)) stop("Chromosome ", chr, " did not produce ", summary_file)
  tab <- fread(summary_file, data.table = TRUE)
  if (nrow(tab) != 1L || !all(c("Chr_num", "ld_window_bp") %in% names(tab)) ||
      is.na(tab$Chr_num[1]) || tab$Chr_num[1] != chr) {
    stop("Unexpected LD-decay summary format for chromosome ", chr)
  }
  summaries[[i]] <- tab
}

windows <- rbindlist(summaries, fill = TRUE, use.names = TRUE)
setorder(windows, Chr_num)
window_file <- file.path(out_dir, "ld_decay_windows_by_chromosome.csv")
fwrite(windows, window_file)

# Do not proceed to LD-block construction with a missing/undefined decay window.
bad_chr <- windows[!is.finite(ld_window_bp) | ld_window_bp <= 0L, Chr_num]
if (length(bad_chr) > 0L) {
  stop("LD decay could not be estimated for chromosome(s): ", paste(bad_chr, collapse = ", "),
       ". Inspect per-chromosome diagnostics in ", out_dir,
       ". The summary was saved to ", window_file)
}

if (all(file.exists(binned_files))) {
  binned <- rbindlist(lapply(binned_files, fread), fill = TRUE, use.names = TRUE)
  setorder(binned, Chr_num, distance_bin)
  fwrite(binned, file.path(out_dir, "ld_decay_binned_by_chromosome.csv"))
}
existing_smooth <- smooth_files[file.exists(smooth_files)]
if (length(existing_smooth) > 0L) {
  smooth <- rbindlist(lapply(existing_smooth, fread), fill = TRUE, use.names = TRUE)
  setorder(smooth, Chr_num, distance_bp)
  fwrite(smooth, file.path(out_dir, "ld_decay_smooth_by_chromosome.csv"))
}

message("Chromosome-specific LD decay calculated: ", window_file)
message("Next: Rscript 6_LD_Block_Concordance.R")
