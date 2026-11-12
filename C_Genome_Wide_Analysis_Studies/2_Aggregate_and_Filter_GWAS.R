#!/usr/bin/env Rscript
# Aggregate GAPIT Filter_GWAS_results from 3-PC GWAS runs.
# IMPORTANT: Filter results are not necessarily BH-FDR q<0.05.
suppressPackageStartupMessages({library(dplyr); library(readr); library(stringr); library(tidyr)})
root <- file.path('outputs', 'GWAS_runs')
out <- file.path('outputs', 'associations')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(root)) stop('No GWAS runs: ', root)
files <- list.files(root, pattern='^GAPIT\\.Association\\.Filter_GWAS_results\\.csv$',
                    full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
if (!length(files)) stop('No GAPIT Filter_GWAS_results.csv files under ', root)
labels <- c(ETH23='results_2023', ETH24='results_2024', ALL='results_all',
            BBCH_ETH23='results_bbch_2023', BBCH_ETH24='results_bbch_2024',
            BBCH_ALL='results_bbch_all')
read_one <- function(f) {
  parent <- basename(dirname(f))
  ds_folder <- basename(dirname(dirname(f)))
  ds <- sub('^results_', '', ds_folder)
  if (!ds %in% names(labels)) stop('Unexpected dataset directory: ', ds_folder)
  prefix <- paste0('GWAS_', ds, '_')
  if (!startsWith(parent, prefix)) stop('Unexpected GWAS directory: ', parent)
  trait_from_dir <- substring(parent, nchar(prefix)+1L)
  x <- read_csv(f, show_col_types=FALSE, name_repair='minimal')
  names(x) <- trimws(names(x))
  if (ncol(x) && (names(x)[1] == '' || grepl('^\\.\\.\\.', names(x)[1]))) x <- x[,-1,drop=FALSE]
  if (!nrow(x)) return(NULL)
  req <- c('SNP','Chr','Pos','P.value')
  if (!all(req %in% names(x))) stop('Required GAPIT fields absent: ', f)
  if (!'MAF' %in% names(x)) x$MAF <- NA_real_
  if (!'traits' %in% names(x)) x$traits <- NA_character_
  x <- x %>% mutate(SNP=as.character(SNP), Chr=as.character(Chr),
                    Pos=as.integer(Pos), P.value=as.numeric(P.value), MAF=as.numeric(MAF),
                    method=str_extract(traits, '^[^.]+'),
                    trait_from_col=str_match(traits, '^[^.]+\\.(.+?)\\(')[,2],
                    pc='3_PC', dataset=unname(labels[ds]),
                    trait=if_else(!is.na(trait_from_col) & nzchar(trait_from_col),
                                  trait_from_col, trait_from_dir),
                    source_file=file.path(ds_folder, parent, basename(f)),
                    SNP_ID=if_else(!is.na(Chr) & !is.na(Pos), paste0(Chr,':',Pos),SNP))
  # Match original aggregate_filter_gwas_3PC.R; do NOT invent BH q-values.
  x <- select(x, -any_of(c('H&B.P.Value','HB.adjusted')))
  x
}
all <- dplyr::bind_rows(lapply(files,read_one))
if (!nrow(all)) stop('All GAPIT Filter tables are empty')
all <- arrange(all, trait, dataset, method, Chr, Pos)
write_csv(all, file.path(out,'all_Filter_GWAS_results_3PC_combined.csv'))
unique <- all %>% group_by(SNP, Chr, Pos, trait) %>%
  summarise(datasets=paste(sort(unique(dataset)),collapse=';'),
            methods=paste(sort(unique(na.omit(method))),collapse=';'),
            min_p=min(P.value,na.rm=TRUE), supporting_tests=n(), .groups='drop')
write_csv(unique,file.path(out,'unique_marker_trait_associations.csv'))
message('Combined Filter results: ',nrow(all),' rows; unique marker-trait pairs: ',nrow(unique))
