# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2023-12-22
#
# Script Description: Prepares the suggestive loci for annotation. 
# Only the top hit, not the whole credible set is used. Positions are liftet to 
# hg19 and the columns renamed to match the pipeline input.
#
# Notes:
# pipeline_name: 12_suggestive_hits_annotation_preparation.R

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
source("helper_scripts/liftOverHg38TOHg19.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

region_data = fread(paste0(path_locus_definition, "locus_definition_suggestive_cleaned.csv"), dec = ",") #TODO check locus definition file
path_output = "09_suggestive_hit_annotation/annotation_input/"

# LIFTING POSITIONS FROM HG38 TO HG19 -------------------------------------

lifted_data = data.table("chr" = region_data$chr, "pos" = region_data$pos, "snps" = region_data$markerID)

lifted_data = liftOVerHg38TOHg19(
  lifted_data,
  try_unlifted_via_ensembl = F,
  return_as_data_table = T,
  path_chainfile = path_chainfile
)

annotation_data = data.table("snp" = lifted_data$snps, "chr_hg19" = lifted_data$chr_hg19, "pos_hg19" = lifted_data$pos_hg19)
annotation_data = merge(annotation_data, region_data[,c("markerID", "aa", "ea", "chrom", "pos", "nWeightedMAF", "region", "phenotype", "betaFEM", "direction", "seFEM", "pFEM", "nWeightedInfoScore", "I2", "totalN")], by.x = "snp", by.y = "markerID", all.x=T, all.y=F )

# CHANGE COLUM NAMES ------------------------------------------------------

setnames(annotation_data, 
         c("chrom", "pos", "nWeightedMAF", "betaFEM", "seFEM", "pFEM", "nWeightedInfoScore", "ea", "aa", "totalN"),
         c("chr_hg38", "pos_hg38", "minor_allele_freq", "effect", "SE", "P-value", "info", "a1", "a2", "n"),
         skip_absent = T)


# OUTPUT ------------------------------------------------------------------
write.table(
  annotation_data,
  paste0(path_output, "input_suggestive_snp_annotation.txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
