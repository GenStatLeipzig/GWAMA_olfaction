# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-09
#
# Script Description: Prepare LDSC command file for genetic correlation
#
#
# Notes:
#
# pipeline_name: R06b_create_ldsc_command.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

p_analysis = "15_ldsr_diseases/"
p_ukbb = "UKBB_phenotypes/"
p_out = paste0(p_analysis, "genetic_correlation/")
f_out = "rg_score_all_diseases"

f_olfaction = "GWASMA_SCORE_all_2024-03-01.sumstats.gz" # can be copied from heritability analysis
files_ukbb = list.files(paste0(p_analysis, p_ukbb), full.names = TRUE)

template = "helper_scripts/ldsr_template_diseases.txt"
script = paste0(p_analysis, "ldsr_command.sh")

directory = p_out
if (!dir.exists(directory)) {
  dir.create(directory)
}


# CREATE COMMAND FILE -----------------------------------------------------

write("", script)
f_olfaction = paste0(p_analysis, f_olfaction)
files = paste0(c(f_olfaction, files_ukbb), collapse = ",")
files = paste0(f_olfaction, ",", files)

ldsr_command = readLines(template)
ldsr_command = str_replace_all(ldsr_command, "FILELIST", files)
ldsr_command = str_replace_all(ldsr_command, "OUT", paste0(p_out, f_out))

write(ldsr_command, script, append = T)
write("\n", script, append = T)

system(paste0("chmod +x ", script))

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
