# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-02-23
#
# Script Description: Extracts a list of genes for each locus that are used for
# candidate gene lookup.
#
#
# Notes:
# pipeline_name: A_02b_extract_gene_list_publication.R
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
source("helper_scripts/liftOverHg38TOHg19.R")
setwd(projectpath)
library("readxl")
library(WriteXLS)

# VARIABLES ---------------------------------------------------------------

annotation_folder = "08_credible_set_analysis/results/" #TODO check folder for annotation results
anno_files = list.files(annotation_folder, pattern = ".xlsx", full.names = T)
out_folder = "08_credible_set_analysis/annotation_summary/" #TODO check output folder
out_file = "gene_list_publication.xlsx"
locus_definition = fread("05_locus_definition/locus_definition_rsID.csv", dec = ",") #TODO check annotation file

# PROCESSING --------------------------------------------------------------

result = foreach (anno_file = anno_files) %do% {
  
  fout = str_split(anno_file, "/", simplify = T)
  fout = fout[length(fout)]
  fout = str_match(fout, ".*_(region.*).xlsx")[2]
  
  # get SNP data
  snp = str_split(fout, "_", n=5, simplify = T)
  snp = snp[length(snp)]
  snp = str_replace_all(snp, "_", ":")
  
  snp_data = locus_definition[markerID == snp,]
  # snp_pos = snp_data[,pos]
  
  # lift snp position to hg19 coordinates to match with annotation
  lifted_data = data.table("chr" = snp_data$chrom, "pos" = snp_data$pos, "snps" = snp_data$markerID)
  
  lifted_data = liftOVerHg38TOHg19(
    lifted_data,
    try_unlifted_via_ensembl = F,
    return_as_data_table = T,
    path_chainfile = path_chainfile
  )
  snp_pos = lifted_data$pos_hg19
  
  if(nrow(snp_data) != 1){
    message("Error when looking up SNP: ", snp, "\n")
    next
  }
  
  
  snp_data_labels = c("region", "phenotype/subgroup", "Top-SNP chromosome", "Top-SNP position (hg38)", "Top-SNP position (hg19)","effect", "se", "pval")
  
  
  snp_data = snp_data[,c("region", "phenotype", "chrom", "pos","betaFEM", "seFEM", "pFEM")]
  snp_data = as.list(snp_data[1])
  snp_data = c(snp_data[1:4], snp_pos, snp_data[5:7])
  snp_data = data.table(info = snp_data_labels, value = snp_data)

  
  # get gene data
  data_genes = as.data.table(read_xlsx(anno_file, sheet = "proximate_genes"))
  gene_list = unique(data_genes$`Abbreviated genename`)
  gene_list = str_sort(gene_list, numeric = T)
  gene_list = c("genes within credible set", gene_list)
  gene_list = na.omit(gene_list)
  
  
  # get the distance of genes to top-snp
  for (row in 1:nrow(data_genes)){
    gene_start = data_genes[row, `start (bp) gene`]
    gene_end = data_genes[row, `end (bp) gene`]
    
    if (is.na(gene_start) | is.na(gene_end)) {
      data_genes[row, top_snp_dist := NA]
    } else if (snp_pos > gene_start & snp_pos < gene_end) {
      data_genes[row, top_snp_dist := 0]
    } else {
      d1 = abs(snp_pos - gene_start)
      d2 = abs(snp_pos - gene_end)
      d = min(d1, d2)
      data_genes[row, top_snp_dist := d]
    }
  }
  
  dist_data = data_genes[`Abbreviated genename` %in% gene_list ,c("Abbreviated genename", "top_snp_dist")]
  dist_data = unique(dist_data)
  setorder(dist_data, top_snp_dist)
  setnames(dist_data, "top_snp_dist", "Distance to index SNP")
  dist_data[,Locus := snp_data$value[[1]]]
  
  dist_data = dist_data[,c("Locus","Abbreviated genename", "Distance to index SNP" )]
  dist_data
  return(dist_data)
  # WriteXLS(
  #   list(snp_data, dist_data, pheno_list),
  #   ExcelFileName = paste0(out_folder, fout, ".xlsx"),
  #   SheetNames = c("snp data", "genes", "phenotypes"),
  #   col.names = T,
  #   row.names = F
  # )
  
}


# SAVE OUTPUT -------------------------------------------------------------
result = rbindlist(result)
setorder(result, Locus, `Distance to index SNP`)
WriteXLS(result,ExcelFileName = paste0(out_folder, out_file))

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
