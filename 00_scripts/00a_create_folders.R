# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-08-16
#
# Script Description: Create folder structure to save script output. If other names 
# are desired, folders need to be created manually and folder names in the 
# scripts need to be adjusted.
#
#
# Notes:
#
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

folder_list = c(
  "01_data_gwcat/",
  "02_MetaGWAS/",
  "03_power_plots/",
  "04_mh_qq_plots_overview_stats/",
  "05_locus_definition/",
  "06_ra_plots/",
  "07a_cojo_loci/input/",
  "07a_cojo_loci/output/",
  "07b_cojo_cond_lit/input/",
  "07b_cojo_cond_lit/output/",
  "07c_sex_conditioning/input/",
  "07c_sex_conditioning/output/",
  "07d_chr11_cojo/input/",
  "07d_chr11_cojo/output/",
  "08_credible_set_analysis/annotation_input/",
  "08_credible_set_analysis/annotation_summary/",
  "08_credible_set_analysis/results/topliste_tabdelim/",
  "09_suggestive_hit_annotation/annotation_input/",
  "09_suggestive_hit_annotation/results/",
  "10_eqtl_coloc/",
  "11_sex_interaction/",
  "12_MR/",
  "12_MR_temp"
)


# PROCESSING --------------------------------------------------------------

for (directory in folder_list){
  if(!dir.exists(directory)){dir.create(directory, recursive = T)}
}

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
