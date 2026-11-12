#!/usr/bin/env Rscript

# GAPIT GWAS using three genomic principal components.

library(GAPIT)
library(dplyr)

fp <- function(...) file.path(...)

load_genotypes <- function() {
  snp_path <- fp("data", "SNP_information.txt")
  gwas_path <- fp("data", "GWAFFs_imputed.txt")
  if (!file.exists(snp_path)) stop("Missing marker metadata: ", snp_path)
  if (!file.exists(gwas_path)) stop("Missing genotype matrix: ", gwas_path)

  myGM <- read.table(snp_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(myGM) < 3L) stop("SNP_information.txt needs SNP, chromosome, position")
  names(myGM)[1:3] <- c("SNP", "Chromosome", "Position")
  myGM$SNP <- as.character(myGM$SNP)
  if (anyDuplicated(myGM$SNP)) stop("Duplicate marker IDs in SNP_information.txt")

  raw <- read.table(gwas_path, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, sep = "", quote = "")
  if (ncol(raw) < 3L) stop("Genotype matrix must have SNP rows and accession columns")
  marker_ids <- as.character(raw[[1]])
  accession_ids <- colnames(raw)[-1]
  if (anyDuplicated(marker_ids)) stop("Duplicate marker IDs in genotype matrix")
  if (anyDuplicated(accession_ids)) stop("Duplicate accession IDs in genotype matrix")
  if (!setequal(marker_ids, myGM$SNP)) {
    stop("SNP_information.txt and genotype matrix have different SNP identifiers")
  }
  dosage <- suppressWarnings(as.matrix(data.frame(lapply(raw[-1], as.numeric),
                                                 check.names = FALSE)) * 2)
  if (anyNA(dosage)) stop("Genotype matrix contains missing/non-numeric values")
  colnames(dosage) <- accession_ids
  rownames(dosage) <- marker_ids
  dosage <- dosage[myGM$SNP, , drop = FALSE]
  gd <- as.data.frame(t(dosage), check.names = FALSE, stringsAsFactors = FALSE)
  gd <- data.frame(taxa = rownames(gd), gd, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(gd) <- NULL
  message("Loaded genotypes: ", nrow(gd), " accessions, ", ncol(gd) - 1L, " markers")
  list(GM = myGM, GD = gd)
}

load_phenotypes <- function(trait, pheno_file) {
  if (!file.exists(pheno_file)) stop("Missing phenotype BLUE file: ", pheno_file)
  df <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
  lower <- tolower(names(df))
  id_index <- match(TRUE, lower %in% c("fe_id", "feid", "genotype.id", "genotype_id", "id", "taxa"))
  blue_index <- which(lower == "blue")[1]
  if (is.na(id_index) || is.na(blue_index)) {
    stop("Phenotype file must have a genotype identifier and BLUE column: ", pheno_file)
  }
  phenotype <- data.frame(taxa = as.character(df[[id_index]]),
                          value = suppressWarnings(as.numeric(df[[blue_index]])),
                          stringsAsFactors = FALSE)
  names(phenotype) <- c("taxa", trait)
  phenotype <- phenotype[!is.na(phenotype[[trait]]) &
                               !is.na(phenotype$taxa) & nzchar(trimws(phenotype$taxa)), , drop = FALSE]
  if (anyDuplicated(phenotype$taxa)) stop("Duplicated phenotype genotype IDs: ", pheno_file)
  if (nrow(phenotype) == 0L) stop("No usable phenotypes in ", pheno_file)
  phenotype
}

run_gwas <- function(dataset, trait, nPCA = 3, phenotype_path) {
  geno_data <- load_genotypes()
  myGM <- geno_data$GM
  myGD <- geno_data$GD
  myY <- load_phenotypes(trait, phenotype_path)

  if (!"taxa" %in% colnames(myGD)) {
    gd_taxa <- rownames(myGD)
  } else {
    gd_taxa <- as.character(myGD$taxa)
  }
  pheno_ids <- as.character(myY$taxa)
  overlap <- intersect(pheno_ids, gd_taxa)
  if (length(overlap) == 0) {
    message("First 20 phenotype IDs: ", paste(head(pheno_ids, 20), collapse=", "))
    message("First 20 genotype taxa: ", paste(head(gd_taxa, 20), collapse=", "))
    stop("No overlap between phenotype IDs and genotype taxa.")
  }

  # Individual GAPIT analysis output (all files under outputs/).
  output_dir <- fp("outputs", "GWAS_runs", paste0("results_", dataset), paste0("GWAS_", dataset, "_", trait))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  original_wd <- getwd()
  setwd(output_dir)
  on.exit({ setwd(original_wd) }, add = TRUE)

  gd_marker_cols <- setdiff(colnames(myGD), "taxa")
  col_classes <- sapply(myGD[gd_marker_cols], function(x) class(x)[1])
  non_numeric_cols <- names(col_classes)[!col_classes %in% c("numeric", "integer")]
  if (length(non_numeric_cols) > 0) {
    message("Found non-numeric genotype columns (will attempt coercion): ", paste(head(non_numeric_cols, 20), collapse = ", "))
    for (cname in non_numeric_cols) {
      original_sample <- head(as.character(myGD[[cname]]), 5)
      myGD[[cname]] <- suppressWarnings(as.numeric(as.character(myGD[[cname]])))
      na_rate <- mean(is.na(myGD[[cname]]))
      message(sprintf("Column %s coerced to numeric; NA rate=%.3f; example values before coercion: %s", cname, na_rate, paste(original_sample, collapse=",")))
      if (na_rate > 0.5) warning(sprintf("High NA rate (%.1f%%) after coercion in column %s", na_rate*100, cname))
    }
  }
  final_classes <- sapply(myGD[gd_marker_cols], function(x) class(x)[1])
  still_non_numeric <- names(final_classes)[!final_classes %in% c("numeric", "integer")]
  if (length(still_non_numeric) > 0) stop("After coercion, some GD columns are still non-numeric: ", paste(still_non_numeric, collapse=", "))

  myGAPIT <- GAPIT(
    Y = myY,
    GD = myGD,
    GM = myGM,
    PCA.total = nPCA,
    model = c("Blink", "FarmCPU", "MLM", "MLMM"),
    Multiple_analysis = TRUE
  )
  # return both GAPIT object and output_dir for caller
  return(list(gapit = myGAPIT, output_dir = output_dir))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript 1.1.GWAS_batch.R <dataset> <trait> <nPCA> <phenotype_CSV_path>")
dataset <- args[1]
trait <- args[2]
nPCA <- as.numeric(args[3])
phenotype_path <- args[4]
if (is.na(nPCA) || nPCA < 0 || nPCA != as.integer(nPCA)) stop("nPCA must be a non-negative integer")

cat("Running GWAS for dataset:", dataset, "and trait:", trait, "with PCA.total=", nPCA, "\n")
result <- run_gwas(dataset, trait, nPCA, phenotype_path)
if (is.list(result) && !is.null(result$output_dir)) {
  cat("GWAS completed. Results saved in ", result$output_dir, "\n")
} else {
  cat("GWAS completed. (output location unknown)\n")
}
