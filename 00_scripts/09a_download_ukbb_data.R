# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-09
#
# Script Description: Download ldsc formatted files for the specified phenotypes
#
#
# Notes:
# pipeline_name: R06a_download_ukbb_data.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)
# VARIABLES ---------------------------------------------------------------

sex <- "male" # all, male, female

if (sex == "all") {
  f_selection = "additional_information/ldsc_selected_phenotypes.xlsx"
  p_out = "UKBB_phenotypes/"
} else if (sex == "female") {
	f_selection = "additional_information/ldsc_selected_phenotypes_females.xlsx"
	p_out = "UKBB_phenotypes_females/"
} else if (sex == "male") {
	f_selection = "additional_information/ldsc_selected_phenotypes_males.xlsx"
	p_out = "UKBB_phenotypes_males/"
}

p_analysis = "15_ldsr_diseases/"



# PROCESSING --------------------------------------------------------------

d = read_xlsx(f_selection) %>% as.data.table()

directory = paste0(p_analysis, p_out)
if (!dir.exists(directory)) {
  dir.create(directory)
}


for (r in 1:nrow(d)) {
  d_row = d[r, ]

  name = d_row$description
  if (str_starts(name, "Non-cancer illness code")) {
    name = str_split_i(name, ": ", 2)
  }
  name = name %>%
    tolower() %>%
    str_replace_all(., " ", "_")

  withr::with_dir(
    paste0(p_analysis, p_out),
    system(d_row$ldsc_sumstat_wget)
  )

  f_ldsc = d_row$ldsc_sumstat_file
  f_new = f_ldsc %>%
    paste0(name, "__", .) %>%
    str_replace_all(., "\\.bgz$", "\\.gz")

  withr::with_dir(
    paste0(p_analysis, p_out),
    file.rename(f_ldsc, f_new)
  )
}

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
