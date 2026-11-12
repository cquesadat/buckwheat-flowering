#!/usr/bin/env Rscript
# Fit spatial models for individual experiments and combine yearly BLUEs by REML.

library(dplyr)

## working directories (relative)
base_dir <- "."
data_dir <- file.path(base_dir, "data")
flower_params_dir <- file.path(base_dir, "flowering_parameters")
out_blue_dir <- file.path(base_dir, "BLUEs_BLUPs")
dir.create(out_blue_dir, recursive = TRUE, showWarnings = FALSE)

## Helpers (copied behaviour from original script)
spats_BLUPs <- function(data_f,col_name){
  require("SpATS")
  data_f$row_f <-   as.factor(data_f$ROW)
  data_f$range_f <-   as.factor(data_f$RANGE)
  data_f$FE_ID <- as.factor(data_f$FE_ID)
  data_f$value <- data_f[[col_name]]
  fit.SpATS_h2 <- SpATS(response = "value", 
                        random = ~  range_f + row_f,
                        spatial = ~ PSANOVA(RANGE, ROW, nseg = c(15,10)),
                        genotype = "FE_ID",
                        genotype.as.random = TRUE,
                        data = data_f,
                        control = list(tolerance = 1e-03,monitor=0)
  )
  h2 <- getHeritability(fit.SpATS_h2)
  BLUPs <- predict(object = fit.SpATS_h2, which = c('FE_ID'))
  BLUPs$h2 <- h2
  BLUPs <- BLUPs[,which(names(BLUPs) %in% c("FE_ID","predicted.values","standard.errors","h2"))]
  return(BLUPs)
}

spats_BLUEs <- function(data_f,col_name){
  require("SpATS")
  data_f$row_f <-   as.factor(data_f$ROW)
  data_f$range_f <-   as.factor(data_f$RANGE)
  data_f$FE_ID <- as.factor(data_f$FE_ID)
  data_f$value <- data_f[[col_name]]
  fit.SpATS_blue <- SpATS(response = "value", 
                          random = ~  range_f + row_f,
                          spatial = ~ PSANOVA(RANGE, ROW, nseg = c(15,10)),
                          genotype = "FE_ID",
                          genotype.as.random = FALSE,
                          data = data_f,
                          control = list(tolerance = 1e-03,monitor=0)
  )
  BLUEs <- predict(object = fit.SpATS_blue, which = c('FE_ID'), return.vcov.matrix = TRUE)
  vcov <- as.matrix(attr(BLUEs, "vcov"))
  weights_vcov <- data.frame(weights_vcov = diag(solve(vcov)))
  BLUEs$weights_vcov <- weights_vcov$weights_vcov
  BLUEs <- BLUEs[,which(names(BLUEs) %in% c("FE_ID","predicted.values","standard.errors","weights_vcov"))]
  return(BLUEs)
}

## fit_REML_extended 
fit_REML_extended <- function(df, parameter = "", m_gen_mat = NULL) {
  require(dplyr)
  print(parameter)
  df <- df %>% mutate(genotype = as.factor(genotype.id), year_site.UID = as.factor(year_site.UID))
  genotype_set <- df %>% ungroup() %>% select(genotype, year_site.UID) %>% unique() %>% group_by(genotype) %>% summarize(n = n()) %>% filter(n > 1) %>% ungroup() %>% select(genotype) %>% unique()
  df <- df %>% filter(genotype %in% genotype_set$genotype)
  df <- df %>% filter(!is.na(BLUE))
  genotype_set <- tidyr::expand_grid(year_site.UID = df$year_site.UID %>% unique(), genotype = df$genotype %>% unique())
  df <- left_join(genotype_set, df, by = c("year_site.UID", "genotype"))
  if (!is.null(m_gen_mat)) {
    genotype_set <- tibble(genotype = as.factor(as.character(rownames(m_gen_mat))))
    df <- inner_join(genotype_set, df, by = "genotype")
    m_gen_mat_ <- m_gen_mat[as.character(df$genotype %>% unique()), as.character(df$genotype %>% unique())]
    m_gen_mat <- m_gen_mat_
    df <- df %>% mutate(genotype = as.factor(genotype))
  }
  df <- droplevels(df)
  df <- df %>% arrange(genotype, year_site.UID)
  weights_ <- df$weights_vcov
  weights_[is.na(df$BLUE)] <- 0
  weights_[is.na(weights_)] <- 0
  n.rep <- mean(with(df, table(genotype)))
  if (!exists("maxiter")) maxiter <<- 1500

  # asreml models
  if (is.null(m_gen_mat)) {
    if(n.rep>2) {
      fit_BLUP <- asreml::asreml(BLUE ~ 1 + year_site.UID, data = df, trace = F, random = ~genotype, residual = ~genotype:diag(year_site.UID), weights = weights_, na.action = list(y = "include", x = "omit"), maxit = maxiter)
      fit_BLUP <- asreml::update.asreml(fit_BLUP)
      if(any(as.data.frame(summary(fit_BLUP)$varcomp)['bound'] == 'B')) {
        fit_BLUP <- asreml::asreml(BLUE ~ 1 + year_site.UID, data = df, trace = F, random = ~genotype, residual = ~genotype:id(year_site.UID), weights = weights_, na.action = list(y = "include", x = "omit"), maxit = maxiter)
        fit_BLUP <- asreml::update.asreml(fit_BLUP)
      }
    } else {
      fit_BLUP <- asreml::asreml(BLUE ~ 1 + year_site.UID, data = df, trace = F, random = ~genotype, residual = ~genotype:id(year_site.UID), weights = weights_, na.action = list(y = "include", x = "omit"), maxit = maxiter)
      fit_BLUP <- asreml::update.asreml(fit_BLUP)
    }
    vc.g <- summary(fit_BLUP)$varcomp['genotype', 'component']
  } else {
    fit_BLUP <- asreml::asreml(BLUE ~ 1 + year_site.UID, data = df, trace = F, random = ~asreml::vm(genotype, m_gen_mat), residual = ~id(genotype):diag(year_site.UID), weights = weights_, na.action = list(y = "include", x = "omit"), maxit = maxiter)
    fit_BLUP <- asreml::update.asreml(fit_BLUP)
    if(any(as.data.frame(summary(fit_BLUP)$varcomp)['bound'] == 'B')) {
      fit_BLUP <- asreml::asreml(BLUE ~ 1 + year_site.UID, data = df, trace = F, random = ~asreml::vm(genotype, m_gen_mat), residual = ~id(genotype):id(year_site.UID), weights = weights_, na.action = list(y = "include", x = "omit"), maxit = maxiter)
      fit_BLUP <- asreml::update.asreml(fit_BLUP)
    }
    vc.g <- summary(fit_BLUP)$varcomp['vm(genotype, m_gen_mat)', 'component']
  }

  if(n.rep>2) {
    fit_BLUE <- asreml::asreml(BLUE ~ 1 + genotype, data = df, trace = F, random = ~year_site.UID, residual = ~genotype:diag(year_site.UID), weights = weights_, na.action = list(y = "include", x = "omit"), maxit = maxiter)
    fit_BLUE <- asreml::update.asreml(fit_BLUE)
  } else {
    fit_BLUE <- asreml::asreml(BLUE ~ 1 + genotype, data = df, trace = F, random = ~year_site.UID, residual = ~genotype:id(year_site.UID), weights = weights_, na.action = list(y = "include", x = "omit"), maxit = maxiter)
    fit_BLUE <- asreml::update.asreml(fit_BLUE)
  }

  vardiff <- function(M){
    n <- length(diag(M))
    Mi <- matrix(diag(M), n,n, byrow = TRUE)
    Mj <- matrix(diag(M), n,n, byrow = FALSE)
    Mdiff <- Mi + Mj -2*M
    return(Mdiff)
  }

  vdBLUP.mat <- predict(fit_BLUP, classify = "genotype", only = "genotype", sed = TRUE)$sed^2
  vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])
  H2Cullis <- 1 - (vdBLUP.avg / (vc.g * 2))
  g.predBLUE <- predict(fit_BLUE, classify = "genotype")
  df_BLUE <- data.frame(g.predBLUE$pvals[, c(1, 2)])
  if (all(is.na(df_BLUE$predicted.value))) {
    g.predBLUE <- predict(fit_BLUE, classify = "genotype", aliased = T)
    df_BLUE <- data.frame(g.predBLUE$pvals[, c(1, 2)])
  }
  names(df_BLUE) <- c("genotype.id", "BLUE")
  df_BLUE <- df_BLUE %>% mutate(genotype.id = as.character(genotype.id))
  g.predBLUP <- predict(fit_BLUP, classify = "genotype")
  df_BLUP <- data.frame(g.predBLUP$pvals[, c(1, 2)])
  names(df_BLUP) <- c("genotype.id", "BLUP")
  df_BLUP <- df_BLUP %>% mutate(genotype.id = as.character(genotype.id))
  df_BLUE_u_BLUP <- inner_join(df_BLUE, df_BLUP, by = "genotype.id")
  H2d.ij.BLUP <- (vardiff(diag(nrow(df_BLUE_u_BLUP)) * vc.g)-vdBLUP.mat)/(vardiff(diag(nrow(df_BLUE_u_BLUP)) * vc.g)) %>% `dimnames<-`(list(df_BLUE_u_BLUP$genotype.id, df_BLUE_u_BLUP$genotype.id))
  H2dBLUP.eq34 <- mean(H2d.ij.BLUP,na.rm = T)
  df_BLUE_u_BLUP$heritability_cullis <- H2Cullis
  df_BLUE_u_BLUP$heritability_H2BLUP <- H2dBLUP.eq34
  df_BLUE_u_BLUP$n.rep <- n.rep
  df_BLUE_u_BLUP$AIC_BLUE <- summary(fit_BLUE)$aic
  df_BLUE_u_BLUP$BIC_BLUE <- summary(fit_BLUE)$bic
  df_BLUE_u_BLUP$AIC_BLUP <- summary(fit_BLUP)$aic
  df_BLUE_u_BLUP$BIC_BLUP <- summary(fit_BLUP)$bic
  return(df_BLUE_u_BLUP)
}


## Processing blocks: mirror original order but use data_dir and out_base
maxiter <- 1500


## Consolidated per-experiment processing using outputs from 1_Flowering_Model_Fitting.R
exps <- c("ETH23", "ETH24", "LfL23", "LfL24", "ETH24L")
for (exp in exps) {
  param_file <- file.path(flower_params_dir, paste0("parameters_", exp, "_clean_das.csv"))
  archived_filenames <- c(
    ETH23 = "parameters_clean_das_2023_ETH.csv",
    ETH24 = "parameters_clean_das_2024_ETH.csv",
    LfL23 = "parameters_clean_das_2023_LfL.csv",
    LfL24 = "parameters_clean_das_2024_LfL.csv",
    ETH24L = "parameters_clean_das_2024_aw.csv"
  )
  if (!file.exists(param_file)) {
    param_file <- file.path(flower_params_dir, archived_filenames[[exp]])
  }
  if (!file.exists(param_file)) {
    stop("Missing cleaned flowering parameters for ", exp)
  }
  design_file <- file.path(data_dir, paste0("Design_", exp, ".csv"))
  if (!file.exists(design_file)) {
    message("Design file missing for ", exp, ": ", design_file)
    next
  }
  design <- read.csv(design_file, stringsAsFactors = FALSE)
  design$Plot <- paste0("Plot", design$ROW, "_", design$RANGE)
  param_data <- read.csv(param_file)
  params <- setdiff(names(param_data), c("Plot", "ID", "Year", "DAS"))
  if (any(grepl("_date$|_doy$", params))) {
    stop("Spatial analysis requires fitted event dates expressed as DAS: ", exp)
  }
  params <- params[vapply(param_data[params], is.numeric, logical(1))]
  out_dir_exp <- file.path(out_blue_dir, exp)
  dir.create(out_dir_exp, recursive = TRUE, showWarnings = FALSE)
  for (param in params) {
    one_param <- param_data[, which(names(param_data) %in% c("Plot","ID",param))]
    one_param <- merge(design, one_param, by = "Plot")

    BLUEs <- spats_BLUEs(one_param, param)
    BLUPs <- spats_BLUPs(one_param, param)
    # SpATS may return estimates in different orders for fixed and random fits.
    key_blue <- as.character(BLUEs[[1]])
    key_blup <- as.character(BLUPs[[1]])
    if (anyDuplicated(key_blue) || anyDuplicated(key_blup) || !setequal(key_blue, key_blup)) {
      stop("BLUE/BLUP genotype identifiers are inconsistent for trait ", param)
    }
    BLUPs <- BLUPs[match(key_blue, key_blup), , drop = FALSE]
    BLUEs_BLUPs <- cbind(BLUEs, BLUPs)
    colnames(BLUEs_BLUPs) <- c("FE_ID","BLUE","std.err","weights_vcov","ID","BLUP","std.err","h2")
    write.csv(BLUEs_BLUPs, file.path(out_dir_exp, paste0(param, "_BLUEs_BLUPs.csv")), row.names = FALSE)
  }
  message("Wrote BLUEs/BLUPs for ", exp, " into ", out_dir_exp)
}


# BBCH parameter CSVs (manual assessments) — write to BLUEs_BLUPs/BBCH_YYYY
bbch_csv_2023 <- file.path(data_dir, "parameter_bbch_2023.csv")
bbch_csv_2024 <- file.path(data_dir, "parameter_bbch_2024.csv")
if (file.exists(bbch_csv_2023)) {
  design <- read.csv(file.path(data_dir, "Design_ETH23.csv"), stringsAsFactors = FALSE)
  design$Plot <- paste0("Plot", design$ROW, "_", design$RANGE)
  param_data <- read.csv(bbch_csv_2023, stringsAsFactors = FALSE)
  if ("genotype_id" %in% names(param_data)) names(param_data)[which(names(param_data)=="genotype_id")] <- "FE_ID"
  if ("genotype.id" %in% names(param_data)) names(param_data)[which(names(param_data)=="genotype.id")] <- "FE_ID"
  out_dir <- file.path(out_blue_dir, "BBCH_2023")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Only BBCH traits: omit CSV row indices and field coordinates.
  metadata <- c("Plot", "FE_ID", "genotype_id", "genotype.id", "row", "range", "ROW", "RANGE", "X", "...1")
  params <- setdiff(names(param_data), metadata)
  params <- params[vapply(param_data[params], is.numeric, logical(1))]
  for (param in params) {
    one_param <- param_data[, which(names(param_data) %in% c("Plot","FE_ID",param))]
    if (!"FE_ID" %in% names(one_param) & "FE_ID" %in% names(param_data)) one_param$FE_ID <- param_data$FE_ID
    colnames(one_param)[which(names(one_param)=="FE_ID")] <- "ID"
    one_param <- merge(design, one_param, by = "Plot")

    BLUEs <- spats_BLUEs(one_param, param)
    BLUPs <- spats_BLUPs(one_param, param)
    # SpATS may return estimates in different orders for fixed and random fits.
    key_blue <- as.character(BLUEs[[1]])
    key_blup <- as.character(BLUPs[[1]])
    if (anyDuplicated(key_blue) || anyDuplicated(key_blup) || !setequal(key_blue, key_blup)) {
      stop("BLUE/BLUP genotype identifiers are inconsistent for trait ", param)
    }
    BLUPs <- BLUPs[match(key_blue, key_blup), , drop = FALSE]
    BLUEs_BLUPs <- cbind(BLUEs, BLUPs)
    colnames(BLUEs_BLUPs) <- c("FE_ID", "BLUE", "std.err", "weights_vcov", "ID", "BLUP", "std.err", "h2")
    write.csv(BLUEs_BLUPs, file.path(out_dir, paste0(param, "_BLUEs_BLUPs.csv")), row.names = FALSE)
  }
}
if (file.exists(bbch_csv_2024)) {
  design <- read.csv(file.path(data_dir, "Design_ETH24.csv"), stringsAsFactors = FALSE)
  design$Plot <- paste0("Plot", design$ROW, "_", design$RANGE)
  param_data <- read.csv(bbch_csv_2024, stringsAsFactors = FALSE)
  if ("genotype_id" %in% names(param_data)) names(param_data)[which(names(param_data)=="genotype_id")] <- "FE_ID"
  if ("genotype.id" %in% names(param_data)) names(param_data)[which(names(param_data)=="genotype.id")] <- "FE_ID"
  out_dir <- file.path(out_blue_dir, "BBCH_2024")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Only BBCH traits: omit CSV row indices and field coordinates.
  metadata <- c("Plot", "FE_ID", "genotype_id", "genotype.id", "row", "range", "ROW", "RANGE", "X", "...1")
  params <- setdiff(names(param_data), metadata)
  params <- params[vapply(param_data[params], is.numeric, logical(1))]
  for (param in params) {
    one_param <- param_data[, which(names(param_data) %in% c("Plot","FE_ID",param))]
    if (!"FE_ID" %in% names(one_param) & "FE_ID" %in% names(param_data)) one_param$FE_ID <- param_data$FE_ID
    colnames(one_param)[which(names(one_param)=="FE_ID")] <- "ID"
    one_param <- merge(design, one_param, by = "Plot")
    BLUEs <- spats_BLUEs(one_param, param)
    BLUPs <- spats_BLUPs(one_param, param)
    # SpATS may return estimates in different orders for fixed and random fits.
    key_blue <- as.character(BLUEs[[1]])
    key_blup <- as.character(BLUPs[[1]])
    if (anyDuplicated(key_blue) || anyDuplicated(key_blup) || !setequal(key_blue, key_blup)) {
      stop("BLUE/BLUP genotype identifiers are inconsistent for trait ", param)
    }
    BLUPs <- BLUPs[match(key_blue, key_blup), , drop = FALSE]
    BLUEs_BLUPs <- cbind(BLUEs, BLUPs)
    colnames(BLUEs_BLUPs) <- c("FE_ID", "BLUE", "std.err", "weights_vcov", "ID", "BLUP", "std.err", "h2")
    write.csv(BLUEs_BLUPs, file.path(out_dir, paste0(param, "_BLUEs_BLUPs.csv")), row.names = FALSE)
  }
}

# Combine across years: follow original approach to build pheno_spats and run fit_REML_extended
combine_and_reml <- function() {
  # read produced BLUEs folders
  read_blue <- function(dir_pattern, env_label) {
    if (!dir.exists(dir_pattern)) return(data.frame())
    files <- list.files(dir_pattern, pattern = "_BLUEs_BLUPs\\.csv$", full.names = TRUE)
    if (length(files)==0) return(data.frame())
    df_list <- lapply(files, function(f) {
      d <- read.csv(f, stringsAsFactors = FALSE)
      if ("FE_ID" %in% names(d)) d$FE_ID <- as.character(d$FE_ID)
      d$trait <- sub("_BLUEs_BLUPs\\.csv$", "", basename(f))
      d$ENV <- env_label
      return(d)
    })
    do.call(rbind, df_list)
  }

  # gather per-year BLUEs
  LfL_2023 <- read_blue(file.path(out_blue_dir, "LfL23"), "LfL_2023")
  LfL_2024 <- read_blue(file.path(out_blue_dir, "LfL24"), "LfL_2024")
  eth_2023 <- read_blue(file.path(out_blue_dir, "ETH23"), "ETH_2023")
  eth_2024 <- read_blue(file.path(out_blue_dir, "ETH24"), "ETH_2024")
  bbch_2023 <- read_blue(file.path(out_blue_dir, "BBCH_2023"), "BBCH_2023")
  bbch_2024 <- read_blue(file.path(out_blue_dir, "BBCH_2024"), "BBCH_2024")

  run_pair_reml <- function(df1, df2, out_subdir) {
    out_dir <- file.path(out_blue_dir, out_subdir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (nrow(df1) == 0 || nrow(df2) == 0) {
      message("Skipping ", out_subdir, " (missing year folder data)")
      return(NULL)
    }
    common_traits <- intersect(unique(df1$trait), unique(df2$trait))
    if (length(common_traits) == 0) {
      message("No common traits to combine for ", out_subdir)
      return(NULL)
    }
    combined <- rbind(df1, df2)
    for (t in common_traits) {
      df_t <- combined[combined$trait == t, , drop = FALSE]
      df_for_reml <- data.frame(
        genotype.id = as.character(df_t$FE_ID),
        BLUE = df_t$BLUE,
        weights_vcov = if ("weights_vcov" %in% names(df_t)) df_t$weights_vcov else NA,
        trait = df_t$trait,
        year_site.UID = df_t$ENV,
        stringsAsFactors = FALSE
      )
      res <- tryCatch({ fit_REML_extended(df_for_reml, parameter = t) }, error = function(e) { message("REML error (", out_subdir, ", ", t, "): ", e$message); NULL })
      if (!is.null(res)) {
        param <- t
        res$trait <- param
        filename <- if (out_subdir == "BBCH_across_years") {
          paste0("BLUPs_BLUEs_", param, ".csv")
        } else {
          paste0(param, "_BLUEs_BLUPs.csv")
        }
        write.csv(res, file.path(out_dir, filename), row.names = FALSE)
      }
    }
  }

  # Explicit across-year combines requested: ETH and BBCH (and LfL optional)
  run_pair_reml(eth_2023, eth_2024, "ETH_across_years")
  run_pair_reml(bbch_2023, bbch_2024, "BBCH_across_years")
  run_pair_reml(LfL_2023, LfL_2024, "LfL_across_years")
}

combine_and_reml()

message("Script finished. Outputs are under: ", out_blue_dir)
