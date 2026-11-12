#!/usr/bin/env Rscript

# Genotype-correlation LD blocks and cross-method GWAS concordance.
# Retains direct SNP overlaps and LD-supported different-position pairs.
# Writes analysis tables; no figures.
#
# LD is computed from the SNP x sample allele-frequency/dosage matrix:
#   r2 = cor(SNP_i, SNP_j)^2

required_packages <- c("tidyverse", "data.table")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with: install.packages(c(",
    paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})



# -------------------------
# Parameters
# -------------------------

r2_threshold <- 0.20
max_block_span_multiplier <- 1.0
script_is_one_level_below_project_root <- FALSE

# -------------------------
# Paths
# -------------------------

get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
  }
  
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  
  getwd()
}

script_dir <- get_script_dir()

project_root <- if (script_is_one_level_below_project_root) {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(script_dir, winslash = "/", mustWork = FALSE)
}

# All generated outputs live under one outputs/ tree.
results_summary <- file.path(project_root, "outputs")

combined_path <- file.path(results_summary, "associations", "all_Filter_GWAS_results_3PC_combined.csv")

snps_info_path <- file.path(project_root, "data", "SNP_information.txt")
geno_path <- file.path(project_root, "data", "GWAFFs_imputed.txt")
ld_windows_path <- file.path(results_summary, "LD_decay", "ld_decay_windows_by_chromosome.csv")

out_dir <- file.path(results_summary, "concordance")
checkpoint_dir <- file.path(out_dir, "chromosome_checkpoints")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

message("GWAS input table: ", normalizePath(combined_path, winslash = "/", mustWork = FALSE))
message("Concordance output directory: ", out_dir)

# -------------------------
# Helpers
# -------------------------

collapse_unique <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[nzchar(x)]
  if (length(x) == 0) return("")
  paste(sort(x), collapse = ";")
}

first_nonmissing <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

safe_r2 <- function(x, y, min_pairs = 5L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_pairs) return(NA_real_)
  
  x <- x[ok]
  y <- y[ok]
  
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  
  as.numeric(cor(x, y)^2)
}

read_genotype_matrix_snp_rows <- function(path) {
  x <- fread(path, data.table = TRUE, check.names = FALSE)
  
  snp_ids <- as.character(x[[1]])
  mat <- as.matrix(x[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  
  rownames(mat) <- snp_ids
  colnames(mat) <- names(x)[-1]
  
  mat
}

get_chr_window <- function(chr_num, ld_windows) {
  w <- ld_windows$ld_window_bp[match(chr_num, ld_windows$Chr_num)]
  
  if (length(w) == 0 || is.na(w) || !is.finite(w) || w <= 0) {
    stop("Missing valid LD-decay window for chromosome: ", chr_num)
  }
  
  as.integer(round(w))
}

pair_category <- function(auto_category, manual_category) {
  case_when(
    auto_category == "flower" & manual_category == "flower" ~ "flower_to_flower",
    auto_category == "abscission" & manual_category == "abscission" ~ "abscission_to_abscission",
    TRUE ~ "mixed"
  )
}

shared_category_records <- function(
    automated_flower,
    automated_abscission,
    manual_flower,
    manual_abscission
) {
  automated_present <- automated_flower == 1L || automated_abscission == 1L
  manual_present <- manual_flower == 1L || manual_abscission == 1L
  
  if (!automated_present || !manual_present) {
    return(character(0))
  }
  
  out <- character(0)
  
  if (automated_flower == 1L && manual_flower == 1L) {
    out <- c(out, "flower_to_flower")
  }
  
  if (automated_abscission == 1L && manual_abscission == 1L) {
    out <- c(out, "abscission_to_abscission")
  }
  
  if (length(out) == 0L) {
    out <- "mixed"
  }
  
  out
}

ld_files_exist <- function(dir_path) {
  all(file.exists(file.path(dir_path, c(
    "snp_to_ld_block.csv",
    "genome_ld_blocks.csv",
    "ld_edges_used.csv"
  ))))
}

copy_ld_files <- function(from_dir, to_dir) {
  files <- c(
    "snp_to_ld_block.csv",
    "genome_ld_blocks.csv",
    "ld_edges_used.csv",
    "ld_windows_used.csv"
  )
  
  for (f in files) {
    src <- file.path(from_dir, f)
    dst <- file.path(to_dir, f)
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
}

# -------------------------
# LD block computation
# -------------------------

make_ld_blocks_one_chr <- function(chr_num, snps_chr, geno_mat, ld_window_bp, max_block_span_bp, r2_threshold) {
  snps_chr <- snps_chr %>% arrange(Pos)
  n <- nrow(snps_chr)
  
  if (n == 0) return(list(snp_blocks = tibble(), edges = tibble(), blocks = tibble()))
  
  parent <- seq_len(n)
  
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  
  union_nodes <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    
    if (ra != rb) {
      parent[rb] <<- ra
    }
  }
  
  edge_rows <- list()
  e <- 1L
  
  pos <- snps_chr$Pos
  ids <- snps_chr$SNP
  
  for (i in seq_len(n)) {
    j <- i + 1L
    
    while (j <= n && (pos[j] - pos[i]) <= ld_window_bp) {
      if (ids[i] %in% rownames(geno_mat) && ids[j] %in% rownames(geno_mat)) {
        r2_val <- safe_r2(as.numeric(geno_mat[ids[i], ]), as.numeric(geno_mat[ids[j], ]))
        
        if (is.finite(r2_val) && r2_val >= r2_threshold) {
          union_nodes(i, j)
          
          edge_rows[[e]] <- tibble(
            Chr_num = chr_num,
            SNP_A = ids[i],
            Pos_A = pos[i],
            SNP_B = ids[j],
            Pos_B = pos[j],
            distance_bp = pos[j] - pos[i],
            raw_r2 = r2_val
          )
          
          e <- e + 1L
        }
      }
      
      j <- j + 1L
    }
  }
  
  roots <- sapply(seq_len(n), find_root)
  comp_id <- match(roots, unique(roots))
  
  snp_tmp <- snps_chr %>%
    mutate(component_id = comp_id)
  
  block_rows <- list()
  snp_rows <- list()
  b <- 1L
  
  for (comp in unique(snp_tmp$component_id)) {
    z <- snp_tmp %>%
      filter(component_id == comp) %>%
      arrange(Pos)
    
    if (nrow(z) == 0) next
    
    # Split very broad connected components to avoid extreme chain effects.
    z <- z %>%
      mutate(chunk = floor((Pos - min(Pos)) / max_block_span_bp) + 1L)
    
    for (chunk_id in unique(z$chunk)) {
      q <- z %>% filter(chunk == chunk_id)
      
      block_id <- paste0("chr", chr_num, "_block_", sprintf("%05d", b))
      
      snp_rows[[b]] <- q %>%
        transmute(
          SNP,
          Chr_num,
          Pos,
          block_id
        )
      
      block_rows[[b]] <- tibble(
        block_id = block_id,
        Chr_num = chr_num,
        block_start = min(q$Pos),
        block_end = max(q$Pos),
        block_length_bp = max(q$Pos) - min(q$Pos) + 1L,
        n_snps = nrow(q),
        member_snps = paste(q$SNP, collapse = ";")
      )
      
      b <- b + 1L
    }
  }
  
  edges <- if (length(edge_rows) == 0) tibble() else bind_rows(edge_rows)
  snp_blocks <- bind_rows(snp_rows)
  blocks <- bind_rows(block_rows)
  
  if (nrow(edges) > 0 && nrow(snp_blocks) > 0) {
    edge_summary <- edges %>%
      left_join(snp_blocks %>% select(SNP_A = SNP, block_A = block_id), by = "SNP_A") %>%
      left_join(snp_blocks %>% select(SNP_B = SNP, block_B = block_id), by = "SNP_B") %>%
      filter(block_A == block_B) %>%
      group_by(block_id = block_A) %>%
      summarise(
        n_ld_edges = n(),
        mean_r2 = mean(raw_r2, na.rm = TRUE),
        min_r2 = min(raw_r2, na.rm = TRUE),
        max_r2 = max(raw_r2, na.rm = TRUE),
        .groups = "drop"
      )
    
    blocks <- blocks %>%
      left_join(edge_summary, by = "block_id")
  } else {
    blocks <- blocks %>%
      mutate(
        n_ld_edges = NA_integer_,
        mean_r2 = NA_real_,
        min_r2 = NA_real_,
        max_r2 = NA_real_
      )
  }
  
  list(snp_blocks = snp_blocks, edges = edges, blocks = blocks)
}

# -------------------------
# Trait groups
# -------------------------

datasets_automated <- c("results_2024", "results_2023", "results_all")
datasets_manual <- c("results_bbch_2024", "results_bbch_2023", "results_bbch_all")

automated_flower_traits <- c(
  "start_flowering_das", "abscission_intercept_amplitude", "area_under_curve",
  "duration_flowering_mid_end", "max_flowering_rate_das", "xmid", "duration_flowering_start_mid",
  "scal", "duration_flowering_start_end", "end_flowering_das"
)

automated_abscission_traits <- c(
  "abscission_slope", "duration_start_flowering_to_end_abscission",
  "abscission_projected_end_das", "abscission_x0",
  "duration_end_flowering_to_end_abscission", "duration_mid_flowering_to_end_abscission"
)

manual_flower_traits <- c(
  "Leaf_1", "Leaf_2", "Leaf_3", "Leaf_4", "Blossoming",
  "First_flower_open", "Flower_open_TINF", "Green_Fruits"
)

manual_abscission_traits <- c(
  "Brown_Fruits", "All_fruits_brown_TINF", "End_Flowering_TINF"
)

# -------------------------
# Read stable input files
# -------------------------

if (!file.exists(ld_windows_path)) stop("Cannot find chromosome-specific LD window file: ", ld_windows_path)

ld_windows_raw <- read_csv(ld_windows_path, show_col_types = FALSE)

if ("ld_window_bp" %in% names(ld_windows_raw)) {
  window_values <- as.numeric(ld_windows_raw$ld_window_bp)
} else if ("background_adjusted_decay_bp" %in% names(ld_windows_raw)) {
  window_values <- as.numeric(ld_windows_raw$background_adjusted_decay_bp)
} else if ("ld_decay_bp" %in% names(ld_windows_raw)) {
  window_values <- as.numeric(ld_windows_raw$ld_decay_bp)
} else {
  stop("LD-window table does not contain ld_window_bp or an equivalent estimate")
}

ld_windows <- tibble(
  Chr_num = as.integer(ld_windows_raw$Chr_num),
  ld_window_bp = window_values
)

if (any(!is.finite(ld_windows$ld_window_bp) | ld_windows$ld_window_bp <= 0)) {
  stop("Invalid chromosome-specific LD-decay windows")
}
write_csv(ld_windows, file.path(out_dir, "ld_windows_used.csv"))

# -------------------------
# Load or compute LD blocks
# -------------------------

if (ld_files_exist(out_dir)) {
  message("Loading existing LD files from output folder: ", out_dir)
  
  snp_to_ld_block <- read_csv(file.path(out_dir, "snp_to_ld_block.csv"), show_col_types = FALSE)
  genome_ld_blocks <- read_csv(file.path(out_dir, "genome_ld_blocks.csv"), show_col_types = FALSE)
  ld_edges <- read_csv(file.path(out_dir, "ld_edges_used.csv"), show_col_types = FALSE)
  
} else {
  message("No complete precomputed LD files found. Computing LD blocks now.")
  
  geno_mat <- read_genotype_matrix_snp_rows(geno_path)
  
  snps_info_raw <- fread(snps_info_path, data.table = FALSE, check.names = FALSE)
  
  snp_col <- grep("^SNP$|marker|Marker", names(snps_info_raw), ignore.case = TRUE, value = TRUE)[1]
  chr_col <- grep("Cromosome|Chromosome|Chr", names(snps_info_raw), ignore.case = TRUE, value = TRUE)[1]
  pos_col <- grep("Position|Pos", names(snps_info_raw), ignore.case = TRUE, value = TRUE)[1]
  
  snps_info <- snps_info_raw %>%
    transmute(
      SNP = as.character(.data[[snp_col]]),
      Chr_num = as.integer(.data[[chr_col]]),
      Pos = as.integer(.data[[pos_col]])
    ) %>%
    filter(SNP %in% rownames(geno_mat)) %>%
    arrange(Chr_num, Pos)
  
  all_snp_blocks <- list()
  all_edges <- list()
  all_blocks <- list()
  
  for (chr in sort(unique(snps_info$Chr_num))) {
    snp_file <- file.path(checkpoint_dir, paste0("snp_to_ld_block_chr", chr, ".csv"))
    edge_file <- file.path(checkpoint_dir, paste0("ld_edges_chr", chr, ".csv"))
    block_file <- file.path(checkpoint_dir, paste0("genome_ld_blocks_chr", chr, ".csv"))
    
    if (file.exists(snp_file) && file.exists(edge_file) && file.exists(block_file)) {
      message("Loading checkpoint for chromosome ", chr, " at ", Sys.time())
      
      all_snp_blocks[[as.character(chr)]] <- read_csv(snp_file, show_col_types = FALSE)
      all_edges[[as.character(chr)]] <- read_csv(edge_file, show_col_types = FALSE)
      all_blocks[[as.character(chr)]] <- read_csv(block_file, show_col_types = FALSE)
      
      next
    }
    
    message("Computing LD blocks for chromosome ", chr, " at ", Sys.time())
    
    chr_window <- get_chr_window(chr, ld_windows)
    max_block_span_bp <- as.integer(round(chr_window * max_block_span_multiplier))
    
    res <- make_ld_blocks_one_chr(
      chr_num = chr,
      snps_chr = snps_info %>% filter(Chr_num == chr),
      geno_mat = geno_mat,
      ld_window_bp = chr_window,
      max_block_span_bp = max_block_span_bp,
      r2_threshold = r2_threshold
    )
    
    write_csv(res$snp_blocks, snp_file)
    write_csv(res$edges, edge_file)
    write_csv(res$blocks, block_file)
    
    message("Saved checkpoint for chromosome ", chr, " at ", Sys.time())
    
    all_snp_blocks[[as.character(chr)]] <- res$snp_blocks
    all_edges[[as.character(chr)]] <- res$edges
    all_blocks[[as.character(chr)]] <- res$blocks
  }
  
  snp_to_ld_block <- bind_rows(all_snp_blocks)
  ld_edges <- bind_rows(all_edges)
  genome_ld_blocks <- bind_rows(all_blocks)
  
  write_csv(snp_to_ld_block, file.path(out_dir, "snp_to_ld_block.csv"))
  write_csv(ld_edges, file.path(out_dir, "ld_edges_used.csv"))
  write_csv(genome_ld_blocks, file.path(out_dir, "genome_ld_blocks.csv"))
}

# -------------------------
# GWAS category points and block comparisons
# -------------------------

# Read as character to avoid parsing chr:pos fields as time.
all_sig <- read_csv(
  combined_path,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

category_points <- all_sig %>%
  filter(pc == "3_PC") %>%
  mutate(
    SNP = as.character(SNP),
    trait = trimws(as.character(trait)),
    P.value = as.numeric(P.value),
    Chr_num = readr::parse_number(as.character(Chr)),
    Pos = as.integer(Pos),
    source = case_when(
      dataset %in% datasets_automated ~ "automated",
      dataset %in% datasets_manual ~ "manual",
      TRUE ~ NA_character_
    ),
    category = case_when(
      source == "automated" & trait %in% automated_flower_traits ~ "flower",
      source == "automated" & trait %in% automated_abscission_traits ~ "abscission",
      source == "manual" & trait %in% manual_flower_traits ~ "flower",
      source == "manual" & trait %in% manual_abscission_traits ~ "abscission",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(source), !is.na(category), !is.na(Chr_num), !is.na(Pos)) %>%
  group_by(SNP, Chr_num, Pos, source, category) %>%
  summarise(
    fisher_stat = sum(-2 * log(pmax(P.value, .Machine$double.xmin)), na.rm = TRUE),
    n_pvalues = sum(!is.na(P.value)),
    fisher_df = 2 * n_pvalues,
    combined_p = pchisq(fisher_stat, df = fisher_df, lower.tail = FALSE),
    logP = -log10(combined_p),
    traits = collapse_unique(trait),
    datasets = collapse_unique(dataset),
    .groups = "drop"
  ) %>%
  mutate(
    x_pos = if_else(source == "automated", 1, 2),
    source_category = paste(source, category, sep = "_")
  ) %>%
  left_join(snp_to_ld_block, by = c("SNP", "Chr_num", "Pos"))

# Pairwise category relations are retained for audit purposes.
# A SNP can have more than one relation if it belongs to multiple trait groups.
direct_hit_pair_relations <- category_points %>%
  filter(source == "automated") %>%
  inner_join(
    category_points %>% filter(source == "manual"),
    by = c("SNP", "Chr_num", "Pos"),
    suffix = c("_auto", "_manual"),
    relationship = "many-to-many"
  ) %>%
  mutate(pair_category = pair_category(category_auto, category_manual)) %>%
  transmute(
    SNP,
    Chr_num,
    Pos,
    auto_category = category_auto,
    manual_category = category_manual,
    pair_category,
    auto_logP = logP_auto,
    manual_logP = logP_manual,
    auto_traits = traits_auto,
    manual_traits = traits_manual,
    block_id = block_id_auto
  )

# Build a four-category marker membership table (one row per marker).
snp_membership <- category_points %>%
  transmute(
    SNP,
    Chr_num,
    Pos,
    set_key = paste(source, category, sep = "_")
  ) %>%
  distinct() %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = set_key,
    values_from = present,
    values_fill = 0L
  )

required_membership_columns <- c(
  "automated_flower",
  "automated_abscission",
  "manual_flower",
  "manual_abscission"
)

for (column_name in required_membership_columns) {
  if (!column_name %in% names(snp_membership)) {
    snp_membership[[column_name]] <- 0L
  }
}

snp_membership <- snp_membership %>%
  select(
    SNP,
    Chr_num,
    Pos,
    all_of(required_membership_columns)
  )

# Stop rather than silently count one SNP at several genomic coordinates.
coordinate_check <- snp_membership %>%
  distinct(SNP, Chr_num, Pos) %>%
  count(SNP, name = "n_coordinates") %>%
  filter(n_coordinates > 1)

if (nrow(coordinate_check) > 0) {
  write_csv(
    coordinate_check,
    file.path(out_dir, "ERROR_SNPs_with_multiple_coordinates.csv")
  )
  stop(
    "At least one SNP occurs at multiple genomic coordinates. See ",
    file.path(out_dir, "ERROR_SNPs_with_multiple_coordinates.csv")
  )
}

direct_metadata <- category_points %>%
  group_by(SNP, Chr_num, Pos) %>%
  summarise(
    auto_logP = if (any(source == "automated")) max(logP[source == "automated"], na.rm = TRUE) else NA_real_,
    manual_logP = if (any(source == "manual")) max(logP[source == "manual"], na.rm = TRUE) else NA_real_,
    auto_traits = collapse_unique(traits[source == "automated"]),
    manual_traits = collapse_unique(traits[source == "manual"]),
    block_id = first_nonmissing(block_id),
    .groups = "drop"
  )

# Derive direct-hit records from the four-category marker membership table.
#
# A SNP produces:
#   - one flower_to_flower record when flowering is present in both sources;
#   - one abscission_to_abscission record when abscission is present in both;
#   - one mixed record only when both sources are represented but no category
#     is shared;
#   - two records when both flowering and abscission are shared.
membership_signature_counts <- snp_membership %>%
  mutate(
    membership_signature = paste0(
      "AF", automated_flower,
      "_AA", automated_abscission,
      "_MF", manual_flower,
      "_MA", manual_abscission
    ),
    automated_any = automated_flower == 1L | automated_abscission == 1L,
    manual_any = manual_flower == 1L | manual_abscission == 1L,
    direct_types = purrr::pmap(
      list(
        automated_flower,
        automated_abscission,
        manual_flower,
        manual_abscission
      ),
      shared_category_records
    ),
    n_direct_records_per_snp = lengths(direct_types),
    direct_types_text = purrr::map_chr(
      direct_types,
      ~paste(.x, collapse = ";")
    )
  ) %>%
  count(
    membership_signature,
    automated_any,
    manual_any,
    n_direct_records_per_snp,
    direct_types_text,
    name = "n_unique_snps"
  ) %>%
  mutate(
    n_direct_records = n_unique_snps * n_direct_records_per_snp
  ) %>%
  arrange(membership_signature)

direct_hits <- snp_membership %>%
  mutate(
    automated_any = automated_flower == 1L | automated_abscission == 1L,
    manual_any = manual_flower == 1L | manual_abscission == 1L,
    has_flower_pair = automated_flower == 1L & manual_flower == 1L,
    has_abscission_pair = automated_abscission == 1L & manual_abscission == 1L,
    has_cross_pair =
      (automated_flower == 1L & manual_abscission == 1L) |
      (automated_abscission == 1L & manual_flower == 1L),
    membership_signature = paste0(
      "AF", automated_flower,
      "_AA", automated_abscission,
      "_MF", manual_flower,
      "_MA", manual_abscission
    ),
    auto_categories = case_when(
      automated_flower == 1L & automated_abscission == 1L ~ "flower;abscission",
      automated_flower == 1L ~ "flower",
      automated_abscission == 1L ~ "abscission",
      TRUE ~ ""
    ),
    manual_categories = case_when(
      manual_flower == 1L & manual_abscission == 1L ~ "flower;abscission",
      manual_flower == 1L ~ "flower",
      manual_abscission == 1L ~ "abscission",
      TRUE ~ ""
    ),
    pair_category = purrr::pmap(
      list(
        automated_flower,
        automated_abscission,
        manual_flower,
        manual_abscission
      ),
      shared_category_records
    )
  ) %>%
  filter(automated_any, manual_any) %>%
  tidyr::unnest_longer(
    pair_category,
    values_to = "pair_category"
  ) %>%
  left_join(
    direct_metadata,
    by = c("SNP", "Chr_num", "Pos")
  ) %>%
  group_by(SNP, Chr_num, Pos) %>%
  mutate(
    n_direct_types_at_snp = n(),
    direct_type_index = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    auto_x = 1,
    manual_x = 2,
    direct_hit_id = paste(SNP, pair_category, sep = "__")
  ) %>%
  select(
    direct_hit_id,
    SNP,
    Chr_num,
    Pos,
    automated_flower,
    automated_abscission,
    manual_flower,
    manual_abscission,
    membership_signature,
    auto_categories,
    manual_categories,
    has_flower_pair,
    has_abscission_pair,
    has_cross_pair,
    pair_category,
    n_direct_types_at_snp,
    direct_type_index,
    auto_x,
    manual_x,
    auto_logP,
    manual_logP,
    auto_traits,
    manual_traits,
    block_id
  )

direct_hits_by_snp <- direct_hits %>%
  group_by(
    SNP,
    Chr_num,
    Pos,
    membership_signature
  ) %>%
  summarise(
    n_direct_hit_records = n(),
    direct_hit_types = collapse_unique(pair_category),
    .groups = "drop"
  )

direct_hit_signature_summary <- direct_hits %>%
  count(
    membership_signature,
    pair_category,
    name = "n_unique_snps"
  ) %>%
  arrange(
    membership_signature,
    pair_category
  )

# Retain the original pairwise direct relations only as an audit table.
# They are not used to count or classify the direct lines in the overview.
direct_hit_pair_relations <- category_points %>%
  filter(source == "automated") %>%
  inner_join(
    category_points %>% filter(source == "manual"),
    by = c("SNP", "Chr_num", "Pos"),
    suffix = c("_auto", "_manual"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    pair_category = pair_category(
      category_auto,
      category_manual
    )
  ) %>%
  transmute(
    SNP,
    Chr_num,
    Pos,
    auto_category = category_auto,
    manual_category = category_manual,
    pair_category,
    auto_logP = logP_auto,
    manual_logP = logP_manual,
    auto_traits = traits_auto,
    manual_traits = traits_manual,
    block_id = block_id_auto
  )

# -------------------------
# LD-supported block records using the same shared-category rule
# -------------------------

# First enumerate all automated/manual pairs assigned to the same LD block.
# A pair is LD-supported only when the two signals occur at DIFFERENT physical
# positions. This is stricter than comparing SNP names because two differently
# formatted SNP identifiers can still describe the same genomic position.
block_pair_relations_all <- category_points %>%
  filter(source == "automated", !is.na(block_id)) %>%
  inner_join(
    category_points %>%
      filter(source == "manual", !is.na(block_id)),
    by = c("Chr_num", "block_id"),
    suffix = c("_auto", "_manual"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    same_physical_position = Pos_auto == Pos_manual,
    pair_category = pair_category(
      category_auto,
      category_manual
    )
  ) %>%
  left_join(
    genome_ld_blocks,
    by = c("block_id", "Chr_num")
  )

# Audit table: these are exact-position overlaps and therefore belong only to
# the direct-hit analysis, never to the LD-supported analysis.
same_position_LD_candidates_excluded <- block_pair_relations_all %>%
  filter(same_physical_position) %>%
  transmute(
    block_id,
    Chr_num,
    block_start,
    block_end,
    pair_category,
    auto_snp = SNP_auto,
    manual_snp = SNP_manual,
    auto_pos = Pos_auto,
    manual_pos = Pos_manual,
    auto_category = category_auto,
    manual_category = category_manual,
    auto_traits = traits_auto,
    manual_traits = traits_manual,
    exclusion_reason = "same physical position: direct hit, not LD-supported"
  )

# Retain only genuinely different-position pairs as LD-supported relations.
block_pair_relations_raw <- block_pair_relations_all %>%
  filter(!same_physical_position) %>%
  transmute(
    block_id,
    Chr_num,
    block_start,
    block_end,
    block_length_bp,
    n_snps,
    n_ld_edges,
    mean_r2,
    min_r2,
    max_r2,
    pair_category,
    auto_snp = SNP_auto,
    manual_snp = SNP_manual,
    auto_pos = Pos_auto,
    manual_pos = Pos_manual,
    auto_category = category_auto,
    manual_category = category_manual,
    auto_logP = logP_auto,
    manual_logP = logP_manual,
    auto_traits = traits_auto,
    manual_traits = traits_manual
  )

block_membership <- category_points %>%
  filter(!is.na(block_id)) %>%
  transmute(
    block_id,
    Chr_num,
    set_key = paste(source, category, sep = "_")
  ) %>%
  distinct() %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = set_key,
    values_from = present,
    values_fill = 0L
  )

for (column_name in required_membership_columns) {
  if (!column_name %in% names(block_membership)) {
    block_membership[[column_name]] <- 0L
  }
}

block_membership <- block_membership %>%
  select(
    block_id,
    Chr_num,
    all_of(required_membership_columns)
  ) %>%
  mutate(
    automated_any = automated_flower == 1L | automated_abscission == 1L,
    manual_any = manual_flower == 1L | manual_abscission == 1L,
    pair_category = purrr::pmap(
      list(
        automated_flower,
        automated_abscission,
        manual_flower,
        manual_abscission
      ),
      shared_category_records
    )
  ) %>%
  filter(automated_any, manual_any) %>%
  tidyr::unnest_longer(
    pair_category,
    values_to = "pair_category"
  ) %>%
  left_join(
    genome_ld_blocks,
    by = c("block_id", "Chr_num")
  )

# Pairwise SNP relations are kept only when they correspond to a block-level
# category record. Cross-category pairings are not labelled mixed when a
# same-category overlap exists in the block.
block_hit_pairs <- block_pair_relations_raw %>%
  inner_join(
    block_membership %>%
      select(
        block_id,
        Chr_num,
        pair_category
      ),
    by = c(
      "block_id",
      "Chr_num",
      "pair_category"
    )
  )

block_source_summary <- category_points %>%
  filter(!is.na(block_id)) %>%
  group_by(
    block_id,
    Chr_num
  ) %>%
  summarise(
    auto_snps = collapse_unique(SNP[source == "automated"]),
    manual_snps = collapse_unique(SNP[source == "manual"]),
    auto_categories = collapse_unique(category[source == "automated"]),
    manual_categories = collapse_unique(category[source == "manual"]),
    auto_traits = collapse_unique(traits[source == "automated"]),
    manual_traits = collapse_unique(traits[source == "manual"]),
    n_auto_snps = n_distinct(SNP[source == "automated"]),
    n_manual_snps = n_distinct(SNP[source == "manual"]),
    .groups = "drop"
  )

block_pair_counts <- block_hit_pairs %>%
  count(
    block_id,
    Chr_num,
    pair_category,
    name = "n_pairwise_comparisons"
  )

block_hits_by_block_all <- block_membership %>%
  left_join(
    block_source_summary,
    by = c("block_id", "Chr_num")
  ) %>%
  left_join(
    block_pair_counts,
    by = c(
      "block_id",
      "Chr_num",
      "pair_category"
    )
  ) %>%
  mutate(
    n_pairwise_comparisons = replace_na(
      n_pairwise_comparisons,
      0L
    )
  ) %>%
  select(
    block_id,
    Chr_num,
    block_start,
    block_end,
    block_length_bp,
    n_snps,
    n_ld_edges,
    mean_r2,
    min_r2,
    max_r2,
    pair_category,
    auto_snps,
    manual_snps,
    auto_categories,
    manual_categories,
    auto_traits,
    manual_traits,
    n_auto_snps,
    n_manual_snps,
    n_pairwise_comparisons
  )

# Any block/category record with zero different-position pairs was created only
# by exact-position/direct overlaps and must not be shown as LD-supported.
block_hits_excluded_no_distinct_positions <- block_hits_by_block_all %>%
  filter(n_pairwise_comparisons == 0L) %>%
  mutate(
    exclusion_reason = "no automated/manual SNP pair at different physical positions"
  )

block_hits_by_block <- block_hits_by_block_all %>%
  filter(n_pairwise_comparisons > 0L) %>%
  arrange(
    Chr_num,
    block_start,
    block_end,
    pair_category
  )

summary_tbl <- tibble(
  r2_threshold = r2_threshold,
  max_block_span_multiplier = max_block_span_multiplier,
  n_ld_blocks = nrow(genome_ld_blocks),
  n_ld_edges = nrow(ld_edges),
  n_category_points = nrow(category_points),
  n_cross_source_upset_snps = n_distinct(direct_hits$SNP),
  n_direct_hit_records = nrow(direct_hits),
  n_SNPs_with_two_direct_hit_records = sum(
    direct_hits_by_snp$n_direct_hit_records == 2L
  ),
  n_direct_hit_pair_relations_audit = nrow(direct_hit_pair_relations),
  n_same_position_LD_pair_candidates_excluded = nrow(
    same_position_LD_candidates_excluded
  ),
  n_LD_block_records_excluded_no_distinct_positions = nrow(
    block_hits_excluded_no_distinct_positions
  ),
  n_block_hit_pairs_distinct_positions = nrow(block_hit_pairs),
  n_block_hit_records_distinct_positions = nrow(block_hits_by_block),
  direct_flower_to_flower = sum(
    direct_hits$pair_category == "flower_to_flower"
  ),
  direct_abscission_to_abscission = sum(
    direct_hits$pair_category == "abscission_to_abscission"
  ),
  direct_mixed = sum(
    direct_hits$pair_category == "mixed"
  ),
  block_flower_to_flower = sum(
    block_hits_by_block$pair_category == "flower_to_flower"
  ),
  block_abscission_to_abscission = sum(
    block_hits_by_block$pair_category == "abscission_to_abscission"
  ),
  block_mixed = sum(
    block_hits_by_block$pair_category == "mixed"
  )
)

# Save all analysis files.
write_csv(
  category_points,
  file.path(out_dir, "category_points_with_blocks.csv")
)
write_csv(
  snp_membership,
  file.path(out_dir, "four_group_SNP_membership.csv")
)
write_csv(
  membership_signature_counts,
  file.path(out_dir, "four_group_membership_signature_counts.csv")
)
write_csv(
  direct_hits,
  file.path(out_dir, "direct_hits_SHARED_CATEGORY_RULE.csv")
)
write_csv(
  direct_hits_by_snp,
  file.path(out_dir, "direct_hits_by_SNP_SHARED_CATEGORY_RULE.csv")
)
write_csv(
  direct_hit_signature_summary,
  file.path(out_dir, "direct_hit_signature_summary_SHARED_CATEGORY_RULE.csv")
)
write_csv(
  direct_hit_pair_relations,
  file.path(out_dir, "direct_hit_pair_relations_AUDIT_ONLY.csv")
)
write_csv(
  block_membership,
  file.path(out_dir, "block_membership_SHARED_CATEGORY_RULE.csv")
)
write_csv(
  same_position_LD_candidates_excluded,
  file.path(out_dir, "same_position_LD_candidates_EXCLUDED.csv")
)
write_csv(
  block_hits_excluded_no_distinct_positions,
  file.path(out_dir, "LD_block_records_EXCLUDED_no_distinct_positions.csv")
)
write_csv(
  block_hit_pairs,
  file.path(out_dir, "block_hit_pairs_DISTINCT_POSITIONS_ONLY.csv")
)
write_csv(
  block_hits_by_block,
  file.path(out_dir, "block_hits_by_block_DISTINCT_POSITIONS_ONLY.csv")
)
write_csv(
  summary_tbl,
  file.path(out_dir, "summary_DISTINCT_POSITION_LD_CORRECTED.csv")
)

# Small text file for meeting notes.
methods_text <- c(
  "LD block concordance pipeline",
  "",
  "1. A fixed genome-wide LD block map was generated from the SNP-by-accession allele-frequency/dosage matrix.",
  "2. Pairwise LD was estimated as r2 = cor(SNP_i, SNP_j)^2.",
  "3. SNP pairs were compared only inside chromosome-specific LD-decay windows.",
  paste0("4. SNP pairs with r2 >= ", r2_threshold, " were retained as LD edges."),
  "5. Connected SNPs were grouped into fixed LD blocks, with broad components split by the chromosome-specific LD window.",
  "6. Concordance reads outputs/associations/all_Filter_GWAS_results_3PC_combined.csv and retains pc == 3_PC.",
  "8. A direct SNP produces one record for each category shared between automated and manual phenotyping.",
  "9. A direct SNP with no shared category but with both sources represented produces one mixed record.",
  "10. A SNP shared for both flowering and abscission produces two direct-hit records at the same genomic position.",
  "11. LD-supported blocks use the same shared-category rule.",
  "12. LD-supported pairs are retained only when automated and manual SNPs have different physical positions; same-position pairs remain direct hits only.",
  "",
  "These are dosage-correlation LD blocks, not formal phased haplotype blocks."
)

writeLines(methods_text, file.path(out_dir, "pipeline_methods_summary.txt"))

print(summary_tbl)
message("Concordance tables saved in: ", out_dir)

