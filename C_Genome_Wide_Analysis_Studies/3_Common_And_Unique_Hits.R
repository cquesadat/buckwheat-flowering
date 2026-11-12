#!/usr/bin/env Rscript
# Produce three CSVs of GWAS hits (common and unique) based on aggregated Filter results.
# Writes exact SNP overlap and source-unique summaries.

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

results_summary <- file.path("outputs", "associations")
out_dir <- file.path("outputs", "concordance")
dir_create(out_dir, recurse = TRUE)
combined_path <- path(results_summary, "all_Filter_GWAS_results_3PC_combined.csv")
if (!file_exists(combined_path)) stop("Combined filter results not found: ", combined_path)

df <- read_csv(combined_path, show_col_types = FALSE)

# Dataset groups (respect the new dataset names; ETH* are automated, BBCH* are measured)
datasets_automated <- c("results_2023", "results_2024", "results_all")
datasets_measured <- c("results_bbch_2023", "results_bbch_2024", "results_bbch_all")

# Trait groups (copied from Overall_density_analysis.r — keep exact names)
group_flower_period_automated<- c("start_flowering_das","abscission_intercept_amplitude","area_under_curve",
                                  "duration_flowering_mid_end", "max_flowering_rate_das","xmid","duration_flowering_start_mid",
                                  "scal","duration_flowering_start_end","end_flowering_das")

group_flower_value_automated <- c("Asym","amplitude_flowering_start_mid","max_flowering_rate_value","amplitude_flowering_mid_end",
                                  "amplitude_flowering_start_end","end_flowering_value")

group_abscission_automated <- c("abscission_slope","duration_start_flowering_to_end_abscission","abscission_projected_end_das","abscission_x0",
                                "duration_end_flowering_to_end_abscission","duration_mid_flowering_to_end_abscission" )

group_early_measured <- c("Cotyledon","Seed")
group_flower_measured <- c("Leaf_1","Leaf_2","Leaf_3","Leaf_4","Brown_Fruits","Blossoming","First_flower_open","Flower_open_TINF","Green_Fruits")
group_abscission_measured <- c("All_fruits_brown_TINF","End_Flowering_TINF")

# Restrict to rows that belong to one of the known datasets
df <- df %>% filter(dataset %in% c(datasets_automated, datasets_measured))

if (nrow(df) == 0) stop("No rows found for specified datasets in combined file.")

# Annotate dataset type and subgroup (traits mapping)
df2 <- df %>% mutate(
  dataset_type = case_when(
    dataset %in% datasets_automated ~ "automated",
    dataset %in% datasets_measured ~ "measured",
    TRUE ~ NA_character_
  ),
  subgroup = case_when(
    trait %in% group_flower_period_automated ~ "flower_period",
    trait %in% group_flower_value_automated ~ "flower_value",
    trait %in% group_abscission_automated ~ "abscission",
    trait %in% group_early_measured ~ "early",
    trait %in% group_flower_measured ~ "flower",
    trait %in% group_abscission_measured ~ "abscission",
    TRUE ~ NA_character_
  )
)

# Aggregate per-SNP summary (min P, methods, datasets, traits)
agg_per_snp <- df2 %>% group_by(SNP, Chr, Pos) %>%
  summarize(
    min_P = min(P.value, na.rm = TRUE),
    datasets = paste(sort(unique(na.omit(dataset))), collapse = ";"),
    dataset_types = paste(sort(unique(na.omit(dataset_type))), collapse = ";"),
    methods = paste(sort(unique(na.omit(method))), collapse = ";"),
    traits = paste(sort(unique(na.omit(trait))), collapse = ";"),
    n_hits = n(),
    .groups = "drop"
  )

# Determine sets per dataset type
automated_snps <- df2 %>% filter(dataset %in% datasets_automated) %>% pull(SNP) %>% unique()
measured_snps <- df2 %>% filter(dataset %in% datasets_measured) %>% pull(SNP) %>% unique()

common_snps <- intersect(automated_snps, measured_snps)
automated_only <- setdiff(automated_snps, common_snps)
measured_only <- setdiff(measured_snps, common_snps)

# Prepare output tables with aggregated info for each set
common_tbl <- agg_per_snp %>% filter(SNP %in% common_snps) %>% arrange(Chr, Pos)
automated_tbl <- agg_per_snp %>% filter(SNP %in% automated_only) %>% arrange(Chr, Pos)
measured_tbl <- agg_per_snp %>% filter(SNP %in% measured_only) %>% arrange(Chr, Pos)

# Output files
out_common <- path(out_dir, "direct_common_hits.csv")
out_auto <- path(out_dir, "automated_unique_hits.csv")
out_meas <- path(out_dir, "measured_unique_hits.csv")

write_csv(common_tbl, out_common)
write_csv(automated_tbl, out_auto)
write_csv(measured_tbl, out_meas)

message("Wrote common hits: ", out_common)
message("Wrote automated-unique hits: ", out_auto)
message("Wrote measured-unique hits: ", out_meas)
