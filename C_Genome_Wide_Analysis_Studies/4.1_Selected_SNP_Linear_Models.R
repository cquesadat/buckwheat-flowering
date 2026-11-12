#!/usr/bin/env Rscript
# Fit BLUE ~ allele_frequency for image-associated SNPs overlapping annotated genes.
# Outputs SNP-trait-gene tables and per-gene model ranking, no plots.
suppressPackageStartupMessages({library(dplyr);library(readr);library(tidyr)})
infile <- 'outputs/associations/annotated_marker_trait_associations.csv'
genofile <- 'data/GWAFFs_imputed.txt'
blue_root <- '../B_Flowering_dynamics_Modelling/BLUEs_BLUPs'
out <- 'outputs/marker_effects'
association_out <- 'outputs/associations/snp_trait_gene_associations.csv'
for (p in c(infile,genofile)) if (!file.exists(p)) stop('Missing: ',p)
if (!dir.exists(blue_root)) stop('Missing Pipeline B BLUEs: ',blue_root)
dir.create(out,recursive=TRUE,showWarnings=FALSE)

ass <- read_csv(infile,show_col_types=FALSE) %>%
  filter(dataset %in% c('results_2023','results_2024','results_all'),in_gene,
         !is.na(gene_ID),!is.na(SNP),!is.na(trait)) %>%
  distinct(gene_ID,SNP,trait,description,hits_CDS)
if (!nrow(ass)) stop('No image-derived SNPs inside genes in annotated association table')

# Read ONLY the selected marker rows; this matrix is SNP x individual, values are
# allele frequencies and MUST NOT be multiplied by 2 for these regression fits.
selected <- unique(ass$SNP)
con <- file(genofile,open='r')
header <- strsplit(trimws(readLines(con,n=1)), '\\s+')[[1]]
if (length(header)<3L) stop('Genotype file has fewer than three columns')
ids <- header[-1]
if (anyDuplicated(ids)) stop('Duplicated genotype individual headers')
snp_values <- list()
repeat {
  lines <- readLines(con,n=2000,warn=FALSE)
  if (!length(lines)) break
  for (line in lines) {
    fields <- strsplit(trimws(line),'\\s+')[[1]]
    if (length(fields)!=length(header)) stop('Inconsistent genotype row width: ',substr(line,1,80))
    if (fields[1] %in% selected) snp_values[[fields[1]]] <- suppressWarnings(as.numeric(fields[-1]))
  }
}
close(con)
missing_snp <- setdiff(selected,names(snp_values))
if (length(missing_snp)) stop('Selected SNPs absent from genotype file: ',paste(head(missing_snp,8),collapse=', '))

get_blue <- local({cache<-new.env(parent=emptyenv()); function(dataset,trait) {
  group<-c('2023'='ETH23','2024'='ETH24','all'='ETH_across_years')[[dataset]]
  f<-file.path(blue_root,group,paste0(trait,'_BLUEs_BLUPs.csv'))
  k<-paste(dataset,trait,sep='::')
  if (exists(k,cache,inherits=FALSE)) return(get(k,cache))
  if (!file.exists(f)) {assign(k,NULL,cache);return(NULL)}
  x<-read_csv(f,show_col_types=FALSE)
  id_col<-intersect(c('FE_ID','genotype.id','genotype_id','ID'),names(x))[1]
  if (is.na(id_col)||!'BLUE'%in%names(x)) stop('Missing genotype ID or BLUE column: ',f)
  x<-x %>% transmute(sample_id=as.character(.data[[id_col]]),BLUE=as.numeric(BLUE))
  assign(k,x,cache);x
}})

fits<-list();i<-0L
for (j in seq_len(nrow(ass))) {
  hit<-ass[j,]
  af<-snp_values[[hit$SNP]]
  snp_df<-data.frame(sample_id=ids,allele_frequency=af,stringsAsFactors=FALSE)
  for (ds in c('2023','2024','all')) {
    blue<-get_blue(ds,hit$trait)
    if (is.null(blue)) next
    merged<-inner_join(blue,snp_df,by='sample_id') %>% filter(is.finite(BLUE),is.finite(allele_frequency))
    n<-nrow(merged); fr<-if (n) mean(merged$allele_frequency) else NA_real_
    maf<-if (is.finite(fr)) min(fr,1-fr) else NA_real_
    o<-list(gene=hit$gene_ID,SNP=hit$SNP,trait=hit$trait,dataset=ds,
            description=hit$description,hits_CDS=hit$hits_CDS,n=n,
            mean_freq=fr,MAF=maf,intercept=NA_real_,slope=NA_real_,
            slope_se=NA_real_,slope_p=NA_real_,r2=NA_real_,adj_r2=NA_real_,
            AIC=NA_real_,BIC=NA_real_,residual_sd=NA_real_)
    if (n>=3 && length(unique(merged$allele_frequency))>=2) {
      fit<-lm(BLUE~allele_frequency,data=merged)
      s<-summary(fit)
      o$intercept<-unname(coef(fit)[1]); o$slope<-unname(coef(fit)[2])
      o$slope_se<-unname(s$coefficients[2,2]);o$slope_p<-unname(s$coefficients[2,4])
      o$r2<-s$r.squared;o$adj_r2<-s$adj.r.squared
      o$AIC<-AIC(fit);o$BIC<-BIC(fit);o$residual_sd<-s$sigma
    }
    i<-i+1L;fits[[i]]<-as_tibble(o)
  }
}
if (!length(fits)) stop('No BLUE files matched the selected SNP-trait pairs')
res<-bind_rows(fits)
write_csv(res,file.path(out,'SNP_fit_summary_all_genes.csv'))

best<-res %>% group_by(gene,trait,dataset) %>%
  arrange(desc(r2),slope_p,.by_group=TRUE) %>% slice_head(n=1) %>% ungroup()
write_csv(best,file.path(out,'SNP_best_per_gene_trait_dataset.csv'))
summary_gene_trait<-best %>% group_by(gene,trait) %>%
  summarise(mean_slope=mean(slope,na.rm=TRUE),mean_slope_se=mean(slope_se,na.rm=TRUE),
            mean_slope_p=mean(slope_p,na.rm=TRUE),mean_r2=mean(r2,na.rm=TRUE),
            total_n=sum(n,na.rm=TRUE),mean_freq=mean(mean_freq,na.rm=TRUE),
            mean_MAF=mean(MAF,na.rm=TRUE),.groups='drop')
write_csv(summary_gene_trait,file.path(out,'SNP_trait_summary_per_gene.csv'))
ranked<-res %>% group_by(gene,SNP) %>%
  summarise(total_n=sum(n,na.rm=TRUE),mean_r2=mean(r2,na.rm=TRUE),
            mean_p=mean(slope_p,na.rm=TRUE),.groups='drop') %>%
  group_by(gene) %>% arrange(desc(mean_r2),mean_p,desc(total_n),.by_group=TRUE) %>%
  slice_head(n=1) %>% ungroup()
write_csv(ranked,file.path(out,'SNP_best_overall_per_gene.csv'))
# One table linking each gene/SNP/trait to 2023/2024/combined linear-fit outputs.
wide<-res %>% select(gene,SNP,trait,description,hits_CDS,dataset,intercept,slope,slope_p,r2) %>%
  mutate(dataset=if_else(dataset=='all','combined',dataset)) %>%
  pivot_wider(names_from=dataset,values_from=c(intercept,slope,slope_p,r2),names_glue='{.value}_{dataset}') %>%
  rename_with(~sub('^intercept_','Intercept_',.),starts_with('intercept_')) %>%
  rename_with(~sub('^slope_p_','p_value_slope_',.),starts_with('slope_p_')) %>%
  rename_with(~sub('^slope_','Slope_',.),starts_with('slope_')) %>%
  rename_with(~sub('^r2_','R2_',.),starts_with('r2_')) %>%
  rename(SNP=SNP,Trait=trait) %>%
  relocate(SNP,Trait,gene,hits_CDS,description)
write_csv(wide,association_out)
message('Saved ',nrow(res),' SNP linear-model fits; ',nrow(ranked),' top gene-level SNPs')
