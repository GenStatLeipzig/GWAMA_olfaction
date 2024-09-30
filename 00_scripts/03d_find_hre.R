# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-04-18
#
# Script Description: Makes a lookup of the candidate gene in references that list
# AREs and EREs. 
#
#
# Notes: Script needs a table with the candidate genes for each locus for lookup (one gene per row)
# pipeline_name: 11d_find_hre.R
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(WriteXLS))


# VARIABLES ---------------------------------------------------------------
#see chromosome X paper for the literature details

gene_annotation_file = "additional_information/gene_table_hre_search.xlsx" # a manually created table cosisting of a list of genes selected for HRE search (e.g. all genes near variants of the CS)

#TODO delete files
ere_file = "additional_information/2004_Bourdeau_EstrogenResponseElements_SupTab1.xlsx" #download here: https://doi.org/10.1210/me.2003-0441
are_file = "additional_information/2016_Wilson_AndrogenResponseElements_SupTables.xlsx" #download here: https://doi.org/10.1038/srep32611
outfile = "11_sex_interaction/hre_annotation.xlsx"

# HRE LOOKUP --------------------------------------------------------------


gene_list = as.data.table(read_xlsx(gene_annotation_file))
are = as.data.table(read_xlsx(are_file, sheet = "SI02_CHIPSEQ_ARE_Annotation"))
are_klf = as.data.table(read_xlsx(are_file, sheet = "SI06_CHIPSEQ_ARE_KLF_MOTIFS"))
ere = as.data.table(read_xlsx(ere_file))

for(row in 1:nrow(gene_list)){
  # look up how many AREs are in the vininity of the genes and what their tier is
  gene = gene_list[row, `Candidate gene`]
  
  if(is.na(gene)){
    next
  }
  
  n_ere = nrow(ere[Human_Gene_Name == gene])
  gene_list[row, EREs := n_ere]
  
  n_are = nrow(are[`Gene Symbol` == gene])
  are_tiers = unique(are[`Gene Symbol` == gene, Tier])
  are_tiers = paste(are_tiers, collapse = ",")
  
  gene_list[row, AREs := n_are]
  gene_list[row, ARE_Tier := are_tiers]
  
  # lookup of AREs that overlap with KLF sites (and check wheather a regulation
  # is confirmed by transcription data)
  n_are_klf = nrow(are_klf[`Gene symbol` == gene])
  is_regulated = are[`Gene Symbol` == gene, `Regulation by AR`]
  is_regulated = any(is_regulated != 0)
  are_tiers = unique(are_klf[`Gene symbol` == gene, `ARE Tier`])
  are_tiers = paste(are_tiers, collapse = ",")
  
  gene_list[row, ARE_KLF := n_are_klf]
  gene_list[row, regulated := is_regulated]
  gene_list[row, ARE_KLF_Tier := are_tiers]
  
}


# SAVING ------------------------------------------------------------------

WriteXLS(gene_list, outfile)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
