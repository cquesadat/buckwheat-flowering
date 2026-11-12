#!/usr/bin/env Rscript
# Annotate the original combined filtered associations with direct gene/CDS hits
# and eggNOG descriptions. This script does not produce figures.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(stringr)})
infile <- 'outputs/associations/all_Filter_GWAS_results_3PC_combined.csv'
gfffile <- 'data/annotated_h2.gff3'
eggfile <- 'data/out.emapper.annotations.tsv'
out <- 'outputs/associations'
for (f in c(infile,gfffile,eggfile)) if (!file.exists(f)) stop('Required input not found: ', f)
dir.create(out,recursive=TRUE,showWarnings=FALSE)

read_eggnog <- function(f) {
  # eggNOG output may include a ## comment preamble and a #query header.
  con <- file(f,'r'); on.exit(close(con))
  first <- readLines(con,n=1000,warn=FALSE)
  pos <- grep('^#?query\\t',first)[1]
  if (is.na(pos)) stop('No query header found in ',f)
  x <- read.delim(f,skip=pos-1L,sep='\t',header=TRUE,check.names=FALSE,
                  comment.char='',quote='',stringsAsFactors=FALSE)
  names(x) <- sub('^#','',names(x))
  if (!'Description' %in% names(x) || !'query' %in% names(x))
    stop('eggNOG annotation must include query and Description columns')
  x %>% transmute(gene_ID=sub('^evm\\.model\\.','evm.TU.',as.character(query)),
                  description=as.character(Description)) %>%
    filter(!is.na(gene_ID),nzchar(gene_ID)) %>% distinct(gene_ID,.keep_all=TRUE)
}

raw <- read_csv(infile,show_col_types=FALSE)
needed <- c('SNP','Chr','Pos','dataset','trait','P.value')
if (!all(needed %in% names(raw))) stop('Combined Filter table lacks: ',paste(setdiff(needed,names(raw)),collapse=', '))
raw <- raw %>% mutate(SNP=as.character(SNP),Chr_num=parse_number(as.character(Chr)),Pos=as.integer(Pos))

gff <- read_tsv(gfffile,comment='#',col_names=FALSE,show_col_types=FALSE,
                col_types=cols(.default='c'))
if (ncol(gff)<9) stop('GFF3 must contain 9 fields')
names(gff)[1:9] <- c('seqid','source','type','start','end','score','strand','phase','attributes')
get_attr <- function(x,key) str_match(x,paste0('(?:^|;)',key,'=([^;]+)'))[,2]
gff <- gff %>% mutate(Chr_num=parse_number(seqid),start=as.integer(start),end=as.integer(end),
                      ID=get_attr(attributes,'ID'),Parent=get_attr(attributes,'Parent'))
genes <- gff %>% filter(type=='gene',!is.na(ID)) %>%
  transmute(Chr_num,start,end,gene_ID=ID,gene_Name=get_attr(attributes,'Name'),
            gene_Annotation=coalesce(get_attr(attributes,'RBH'),get_attr(attributes,'BBH')),
            gene_strand=strand)
transcripts <- gff %>% filter(type %in% c('mRNA','transcript'),!is.na(ID)) %>%
  transmute(transcript_id=ID,gene_ID=Parent)
cds <- gff %>% filter(type=='CDS',!is.na(Parent)) %>%
  select(Chr_num,start,end,Parent) %>% separate_rows(Parent,sep=',') %>%
  left_join(transcripts,by=c('Parent'='transcript_id')) %>%
  mutate(gene_ID=coalesce(gene_ID,Parent)) %>% select(Chr_num,start,end,gene_ID)
markers <- raw %>% distinct(SNP,Chr_num,Pos) %>% filter(!is.na(Chr_num),!is.na(Pos))
within <- markers %>% inner_join(genes,by='Chr_num',relationship='many-to-many') %>%
  filter(Pos>=start,Pos<=end) %>% transmute(SNP,Chr_num,Pos,gene_ID,gene_Name,gene_Annotation,
            gene_start=start,gene_end=end,gene_strand) %>% distinct()
in_cds <- markers %>% inner_join(cds,by='Chr_num',relationship='many-to-many') %>%
  filter(Pos>=start,Pos<=end) %>% select(SNP,Chr_num,Pos,gene_ID) %>% distinct() %>% mutate(hits_CDS=TRUE)
anno <- read_eggnog(eggfile)
gene_map <- within %>% left_join(in_cds,by=c('SNP','Chr_num','Pos','gene_ID')) %>%
  mutate(hits_CDS=coalesce(hits_CDS,FALSE),in_gene=TRUE) %>%
  left_join(anno,by='gene_ID') %>% arrange(Chr_num,Pos,gene_ID)
write_csv(gene_map,file.path(out,'marker_gene_overlaps.csv'))
annotated <- raw %>% left_join(gene_map,by=c('SNP','Chr_num','Pos'),relationship='many-to-many') %>%
  mutate(in_gene=coalesce(in_gene,FALSE),hits_CDS=coalesce(hits_CDS,FALSE))
write_csv(annotated,file.path(out,'annotated_marker_trait_associations.csv'))

# One row per gene, image-derived hits, matching the original analysis fields.
# Origin is the phenotypic trait category, NOT a year or experimental site.
period <- c('start_flowering_das','abscission_intercept_amplitude','area_under_curve',
            'duration_flowering_mid_end','max_flowering_rate_das','xmid',
            'duration_flowering_start_mid','scal','duration_flowering_start_end','end_flowering_das')
value <- c('Asym','amplitude_flowering_start_mid','max_flowering_rate_value',
           'amplitude_flowering_mid_end','amplitude_flowering_start_end','end_flowering_value')
abs <- c('abscission_slope','duration_start_flowering_to_end_abscission',
         'abscission_projected_end_das','abscission_x0',
         'duration_end_flowering_to_end_abscission','duration_mid_flowering_to_end_abscission')
image_associations <- annotated %>% filter(dataset %in% c('results_2023','results_2024','results_all')) %>%
  mutate(subgroup=case_when(trait %in% period~'flower_period',
                            trait %in% value~'flower_value',
                            trait %in% abs~'abscission',TRUE~NA_character_),
         logP=-log10(pmax(as.numeric(P.value),.Machine$double.xmin))) %>%
  filter(!is.na(subgroup))
image_gene_rows <- image_associations %>% filter(in_gene) %>%
  group_by(SNP,Chr_num,Pos,gene_ID,gene_Name,gene_Annotation,gene_start,gene_end,gene_strand,
           hits_CDS) %>%
  summarise(Origin=paste(sort(unique(paste0('automated_',subgroup))),collapse=';'),
            traits=paste(sort(unique(trait)),collapse=';'),
            logP_sum=sum(logP,na.rm=TRUE),.groups='drop')
auto <- image_gene_rows %>% group_by(gene_ID,gene_Name,gene_Annotation,
                 gene_start,gene_end,gene_strand) %>%
  summarise(SNP=paste(sort(unique(SNP)),collapse=';'),
            Origin=paste(sort(unique(Origin)),collapse=';'),
            traits=paste(sort(unique(traits)),collapse=';'),
            logP_sum_collapsed=paste(unique(as.character(logP_sum)),collapse=';'),
            logP_sum_sum=sum(logP_sum,na.rm=TRUE),
            hit_CDS=any(hits_CDS),.groups='drop')
write_csv(auto,file.path(out,'automated_unique_hits.csv'))

# Original analysis statistics for each automated phenotype category.
stat_rows <- lapply(c('flower_period','flower_value','abscission','all'),function(g) {
  subset <- if (g=='all') image_associations else
    image_associations %>% filter(subgroup==g)
  snps <- subset %>% distinct(SNP)
  mapped <- gene_map %>% filter(SNP %in% snps$SNP)
  tibble(subgroup=g,initial_unique_snps=nrow(snps),
         snps_in_gene=n_distinct(mapped$SNP),unique_genes_hit=n_distinct(mapped$gene_ID),
         genes_with_cds_hit=n_distinct(mapped$gene_ID[mapped$hits_CDS]))
})
write_csv(bind_rows(stat_rows),file.path(out,'automated_snps_summary_stats.csv'))

message('Mapped ',n_distinct(gene_map$SNP),' markers to genes; ',nrow(auto),' image-associated genes')
