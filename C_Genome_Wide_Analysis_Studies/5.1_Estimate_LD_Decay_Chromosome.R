#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  getwd()
}

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  p <- list(
    chr = NA_integer_,
    max_distance_bp = 2000000L,
    bin_width_bp = 5000L,
    min_pairs_per_bin = 100L,
    loess_span = 0.45,
    plot_grid_bp = 1000L,
    initial_start_bp = 5000L,
    initial_end_bp = 50000L,
    background_start_bp = 500000L,
    background_end_bp = 2000000L,
    consecutive_points = 3L,
    min_window_bp = 50000L
  )
  i <- 1L
  while (i <= length(a)) {
    k <- a[i]
    v <- if (i < length(a)) a[i + 1L] else NA_character_
    if (k == "--chr") p$chr <- as.integer(v)
    if (k == "--max-distance-bp") p$max_distance_bp <- as.integer(v)
    if (k == "--bin-width-bp") p$bin_width_bp <- as.integer(v)
    if (k == "--min-pairs-per-bin") p$min_pairs_per_bin <- as.integer(v)
    if (k == "--loess-span") p$loess_span <- as.numeric(v)
    if (k == "--plot-grid-bp") p$plot_grid_bp <- as.integer(v)
    if (k == "--initial-start-bp") p$initial_start_bp <- as.integer(v)
    if (k == "--initial-end-bp") p$initial_end_bp <- as.integer(v)
    if (k == "--background-start-bp") p$background_start_bp <- as.integer(v)
    if (k == "--background-end-bp") p$background_end_bp <- as.integer(v)
    if (k == "--consecutive-points") p$consecutive_points <- as.integer(v)
    if (k == "--min-window-bp") p$min_window_bp <- as.integer(v)
    i <- i + 2L
  }
  p
}

args <- parse_args()
if (is.na(args$chr)) stop("Missing chromosome. Use: --chr 1")
if (args$background_end_bp > args$max_distance_bp) args$background_end_bp <- args$max_distance_bp

message("Chromosome: ", args$chr)
message("LD measure: raw r2")
message("Mode: all local pairs + background-adjusted half-decay")
message("max_distance_bp: ", args$max_distance_bp)
message("bin_width_bp: ", args$bin_width_bp)
message("initial range: ", args$initial_start_bp, "-", args$initial_end_bp, " bp")
message("background range: ", args$background_start_bp, "-", args$background_end_bp, " bp")
message("minimum downstream window: ", args$min_window_bp, " bp")

script_dir <- get_script_dir()
project_root <- script_dir
results_summary <- file.path(project_root, "outputs")

snps_info_path <- file.path(project_root, "data", "SNP_information.txt")
geno_path <- file.path(project_root, "data", "GWAFFs_imputed.txt")

out_root <- file.path(results_summary, "LD_decay")
out_chr <- file.path(out_root, paste0("chr", args$chr))
dir.create(out_chr, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) stop("Cannot find ", label, ": ", path)
}
stop_if_missing(snps_info_path, "SNP information file")
stop_if_missing(geno_path, "genotype/allele-frequency file")

read_genotype_matrix_snp_rows <- function(path) {
  message("Reading genotype/allele-frequency file...")
  dt <- fread(path, check.names = FALSE, data.table = TRUE)
  if (ncol(dt) < 3L) stop("Genotype file must have one SNP ID column plus individual columns")
  snp_ids <- as.character(dt[[1]])
  if (anyDuplicated(snp_ids)) {
    dup <- unique(snp_ids[duplicated(snp_ids)])
    stop("Duplicated SNP IDs. First duplicated IDs: ", paste(head(dup, 10), collapse = ", "))
  }
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- snp_ids
  colnames(mat) <- names(dt)[-1]
  message("  Loaded ", nrow(mat), " SNPs x ", ncol(mat), " individuals")
  if (anyNA(mat)) {
    message("  Missing values detected. Imputing missing values by SNP mean.")
    snp_means <- rowMeans(mat, na.rm = TRUE)
    bad <- !is.finite(snp_means)
    if (any(bad)) stop("Some SNPs are entirely missing. Remove them before running.")
    idx <- which(is.na(mat), arr.ind = TRUE)
    mat[idx] <- snp_means[idx[, 1]]
  }
  mat
}

normalize_rows <- function(mat) {
  mat <- mat - rowMeans(mat)
  row_norms <- sqrt(rowSums(mat^2))
  keep <- is.finite(row_norms) & row_norms > 0
  mat <- mat[keep, , drop = FALSE]
  row_norms <- row_norms[keep]
  mat <- mat / row_norms
  list(mat = mat, keep = keep)
}

count_pairs_per_bin <- function(pos, max_distance_bp, bin_width_bp, n_bins) {
  n <- length(pos)
  counts <- integer(n_bins)
  total_pairs <- 0
  for (i in seq_len(n - 1L)) {
    j_end <- findInterval(pos[i] + max_distance_bp, pos)
    if (j_end <= i) next
    j <- (i + 1L):j_end
    d <- pos[j] - pos[i]
    bin_idx <- floor(d / bin_width_bp) + 1L
    bin_idx <- bin_idx[bin_idx >= 1L & bin_idx <= n_bins]
    if (length(bin_idx) > 0) {
      counts <- counts + tabulate(bin_idx, nbins = n_bins)
      total_pairs <- total_pairs + length(bin_idx)
    }
  }
  list(counts = counts, total_pairs = total_pairs)
}

fill_r2_bins <- function(norm_mat, pos, max_distance_bp, bin_width_bp, counts) {
  n_bins <- length(counts)
  r2_bins <- vector("list", n_bins)
  for (b in seq_len(n_bins)) r2_bins[[b]] <- if (counts[b] > 0L) numeric(counts[b]) else numeric(0)
  ptr <- integer(n_bins)
  n <- length(pos)
  for (i in seq_len(n - 1L)) {
    j_end <- findInterval(pos[i] + max_distance_bp, pos)
    if (j_end <= i) next
    j <- (i + 1L):j_end
    d <- pos[j] - pos[i]
    bin_idx <- floor(d / bin_width_bp) + 1L
    keep <- bin_idx >= 1L & bin_idx <= n_bins
    if (!any(keep)) next
    j <- j[keep]
    bin_idx <- bin_idx[keep]
    r <- as.numeric(norm_mat[j, , drop = FALSE] %*% norm_mat[i, ])
    r2 <- r^2
    valid <- is.finite(r2)
    if (!any(valid)) next
    bin_idx <- bin_idx[valid]
    r2 <- r2[valid]
    split_vals <- split(r2, bin_idx)
    for (b_chr in names(split_vals)) {
      b <- as.integer(b_chr)
      vals <- split_vals[[b_chr]]
      n_new <- length(vals)
      start <- ptr[b] + 1L
      end <- ptr[b] + n_new
      r2_bins[[b]][start:end] <- vals
      ptr[b] <- end
    }
    if (i %% 500L == 0L) message("  Processed SNP ", i, " / ", n)
  }
  for (b in seq_len(n_bins)) {
    if (ptr[b] < length(r2_bins[[b]])) {
      r2_bins[[b]] <- if (ptr[b] == 0L) numeric(0) else r2_bins[[b]][seq_len(ptr[b])]
    }
  }
  r2_bins
}

make_binned_table <- function(r2_bins, chr_num, chromosome_label, bin_width_bp) {
  rows <- vector("list", length(r2_bins))
  for (b in seq_along(r2_bins)) {
    vals <- r2_bins[[b]]
    vals <- vals[is.finite(vals)]
    bin_start <- (b - 1L) * bin_width_bp
    bin_end <- b * bin_width_bp
    bin_mid <- bin_start + bin_width_bp / 2
    rows[[b]] <- data.table(
      Chr_num = chr_num,
      chromosome = chromosome_label,
      distance_bin = bin_mid,
      bin_start_bp = bin_start,
      bin_end_bp = bin_end,
      n_pairs = length(vals),
      ld_mean = if (length(vals) > 0L) mean(vals) else NA_real_,
      ld_median = if (length(vals) > 0L) median(vals) else NA_real_
    )
  }
  rbindlist(rows)
}

first_consecutive_crossing <- function(dt, threshold, start_bp, n_consecutive) {
  y <- copy(dt[distance_bp >= start_bp & is.finite(ld_smooth)])
  if (nrow(y) == 0) return(list(distance = NA_real_, status = "no_points"))
  y[, below := ld_smooth <= threshold]
  if (n_consecutive <= 1L) {
    val <- y[below == TRUE, distance_bp[1]]
    status <- "single_point_crossing"
  } else {
    run <- rep(TRUE, nrow(y))
    for (lag_i in 0:(n_consecutive - 1L)) {
      run <- run & shift(y$below, lag_i, fill = FALSE)
    }
    y[, below_run := run]
    val <- y[below_run == TRUE, distance_bp[1]]
    status <- paste0(n_consecutive, "_point_crossing")
  }
  if (length(val) == 0L || !is.finite(val)) {
    val <- y[below == TRUE, distance_bp[1]]
    status <- "single_point_crossing_fallback"
  }
  if (length(val) == 0L || !is.finite(val)) {
    val <- max(y$distance_bp, na.rm = TRUE)
    status <- "threshold_not_reached"
  }
  list(distance = val, status = status)
}

fit_and_summarise <- function(binned, chr_num, chromosome_label, n_markers, total_pairs) {
  fit_df <- binned[n_pairs >= args$min_pairs_per_bin & is.finite(ld_median)]
  base <- data.table(
    Chr_num = chr_num,
    chromosome = chromosome_label,
    n_markers = n_markers,
    n_pairs = total_pairs,
    ld_measure = "raw_r2",
    max_distance_bp = args$max_distance_bp,
    bin_width_bp = args$bin_width_bp,
    initial_start_bp = args$initial_start_bp,
    initial_end_bp = args$initial_end_bp,
    background_start_bp = args$background_start_bp,
    background_end_bp = args$background_end_bp,
    min_window_bp = args$min_window_bp
  )
  if (nrow(fit_df) < 8L) {
    return(list(
      summary = cbind(base, data.table(
        naive_max_ld = NA_real_,
        naive_halfmax_threshold = NA_real_,
        naive_decay_bp = NA_real_,
        initial_ld = NA_real_,
        background_ld = NA_real_,
        background_adjusted_threshold = NA_real_,
        background_adjusted_decay_bp = NA_real_,
        ld_decay_bp = NA_real_,
        ld_decay_kb = NA_real_,
        ld_window_bp = NA_real_,
        ld_window_kb = NA_real_,
        fit_status = "too_few_bins"
      )),
      smooth = data.table()
    ))
  }

  fit <- loess(
    ld_median ~ distance_bin,
    data = fit_df,
    weights = sqrt(fit_df$n_pairs),
    span = args$loess_span,
    control = loess.control(surface = "direct")
  )
  smooth <- data.table(distance_bp = seq(0, args$max_distance_bp, by = args$plot_grid_bp))
  smooth[, `:=`(
    Chr_num = chr_num,
    chromosome = chromosome_label,
    ld_smooth = as.numeric(predict(fit, newdata = data.frame(distance_bin = distance_bp)))
  )]
  smooth <- smooth[is.finite(ld_smooth)]
  setorder(smooth, distance_bp)

  naive_max_ld <- max(smooth$ld_smooth, na.rm = TRUE)
  naive_half <- naive_max_ld / 2
  naive_peak <- smooth[which.max(ld_smooth), distance_bp]
  naive_cross <- first_consecutive_crossing(
    dt = smooth,
    threshold = naive_half,
    start_bp = naive_peak,
    n_consecutive = args$consecutive_points
  )

  initial_ld <- fit_df[
    bin_start_bp >= args$initial_start_bp & bin_end_bp <= args$initial_end_bp,
    median(ld_median, na.rm = TRUE)
  ]
  if (!is.finite(initial_ld)) {
    initial_ld <- fit_df[
      distance_bin >= args$initial_start_bp & distance_bin <= args$initial_end_bp,
      median(ld_median, na.rm = TRUE)
    ]
  }

  background_ld <- fit_df[
    bin_start_bp >= args$background_start_bp & bin_end_bp <= args$background_end_bp,
    median(ld_median, na.rm = TRUE)
  ]
  if (!is.finite(background_ld)) {
    q75 <- quantile(fit_df$distance_bin, 0.75, na.rm = TRUE)
    background_ld <- fit_df[distance_bin >= q75, median(ld_median, na.rm = TRUE)]
  }

  adjusted_threshold <- background_ld + 0.5 * (initial_ld - background_ld)
  adjusted_cross <- first_consecutive_crossing(
    dt = smooth,
    threshold = adjusted_threshold,
    start_bp = args$initial_start_bp,
    n_consecutive = args$consecutive_points
  )

  ld_decay_bp <- adjusted_cross$distance
  ld_window_bp <- ifelse(is.finite(ld_decay_bp), max(ld_decay_bp, args$min_window_bp), NA_real_)

  summary <- cbind(base, data.table(
    naive_max_ld = naive_max_ld,
    naive_halfmax_threshold = naive_half,
    naive_decay_bp = naive_cross$distance,
    naive_fit_status = naive_cross$status,
    initial_ld = initial_ld,
    background_ld = background_ld,
    background_adjusted_threshold = adjusted_threshold,
    background_adjusted_decay_bp = adjusted_cross$distance,
    background_adjusted_fit_status = adjusted_cross$status,
    ld_decay_bp = ld_decay_bp,
    ld_decay_kb = ld_decay_bp / 1000,
    ld_window_bp = ld_window_bp,
    ld_window_kb = ld_window_bp / 1000,
    fit_status = "ok"
  ))
  list(summary = summary, smooth = smooth)
}

message("Reading SNP information...")
snps_info <- fread(snps_info_path, data.table = TRUE)
snp_col <- grep("^SNP$|marker|Marker", names(snps_info), ignore.case = TRUE, value = TRUE)[1]
chr_col <- grep("Cromosome|Chromosome|Chr", names(snps_info), ignore.case = TRUE, value = TRUE)[1]
pos_col <- grep("Position|Pos", names(snps_info), ignore.case = TRUE, value = TRUE)[1]
if (is.na(snp_col) || is.na(chr_col) || is.na(pos_col)) {
  stop("SNP_information.txt must contain SNP, chromosome and position columns")
}
snps_info <- snps_info[, .(
  SNP = as.character(get(snp_col)),
  Chr_num = as.integer(get(chr_col)),
  Pos = as.integer(get(pos_col))
)]
snps_info <- snps_info[!is.na(SNP) & !is.na(Chr_num) & !is.na(Pos)]

# Chromosome labels use SNP_information.txt; assembly lengths are not required for LD decay.
chr_label <- paste0("chr", args$chr)

snp_chr <- snps_info[Chr_num == args$chr]
setorder(snp_chr, Pos)
if (nrow(snp_chr) < 2L) stop("Chromosome ", args$chr, " has fewer than 2 SNPs in SNP_information.txt")

geno_mat <- read_genotype_matrix_snp_rows(geno_path)
snp_chr <- snp_chr[SNP %in% rownames(geno_mat)]
if (nrow(snp_chr) < 2L) stop("Chromosome ", args$chr, " has fewer than 2 SNPs present in genotype file")
setorder(snp_chr, Pos)
geno_chr <- geno_mat[snp_chr$SNP, , drop = FALSE]
message("  Chr ", chr_label, ": ", nrow(geno_chr), " SNPs x ", ncol(geno_chr), " individuals")

message("Normalizing SNP rows...")
norm_res <- normalize_rows(geno_chr)
geno_norm <- norm_res$mat
snp_chr <- snp_chr[norm_res$keep]
if (nrow(geno_norm) < 2L) stop("Fewer than two non-constant SNPs remain on chromosome ", args$chr)
message("  Non-constant SNPs kept: ", nrow(geno_norm))

n_bins <- ceiling(args$max_distance_bp / args$bin_width_bp)

message("Counting all SNP pairs within ", args$max_distance_bp, " bp...")
count_res <- count_pairs_per_bin(
  pos = snp_chr$Pos,
  max_distance_bp = args$max_distance_bp,
  bin_width_bp = args$bin_width_bp,
  n_bins = n_bins
)
message("Total pairs to compute: ", count_res$total_pairs)
if (count_res$total_pairs == 0) stop("No SNP pairs found within max_distance_bp")

message("Computing raw r2 for all local pairs...")
r2_bins <- fill_r2_bins(
  norm_mat = geno_norm,
  pos = snp_chr$Pos,
  max_distance_bp = args$max_distance_bp,
  bin_width_bp = args$bin_width_bp,
  counts = count_res$counts
)

binned <- make_binned_table(
  r2_bins = r2_bins,
  chr_num = args$chr,
  chromosome_label = chr_label,
  bin_width_bp = args$bin_width_bp
)

fit_res <- fit_and_summarise(
  binned = binned,
  chr_num = args$chr,
  chromosome_label = chr_label,
  n_markers = nrow(geno_norm),
  total_pairs = sum(binned$n_pairs, na.rm = TRUE)
)

summary_file <- file.path(out_chr, paste0("ld_decay_window_chr", args$chr, ".csv"))
binned_file <- file.path(out_chr, paste0("ld_decay_binned_chr", args$chr, ".csv"))
smooth_file <- file.path(out_chr, paste0("ld_decay_smooth_chr", args$chr, ".csv"))
fwrite(fit_res$summary, summary_file)
fwrite(binned, binned_file)
if (nrow(fit_res$smooth) > 0) fwrite(fit_res$smooth, smooth_file)
message("LD decay summary: ", summary_file)
