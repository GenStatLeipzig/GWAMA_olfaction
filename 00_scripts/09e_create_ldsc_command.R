# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-09
#
# Script Description: Preapare LDSC command file for genetic correlation
#
#
# Notes:
# Script file has to be executed manually.
#
# pipeline_name: R12b_create_ldsc_command.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

p_analysis = "18_ldsr_coffee/"
p_ukbb = "sum_stats_formatted/"
p_out = paste0(p_analysis, "genetic_correlation/")
f_out_template = "rg_PHENO_beverage_consumption"

fl_olfaction = c(
	"GWASMA_coffee_all_2024-03-01.sumstats.gz" # needs to be manualls copied from heritability calculation
) 
files_ukbb = list.files(paste0(p_analysis, p_ukbb), pattern = "\\.gz$", full.names = TRUE)

template = "helper_scripts/ldsr_template.txt"
script = paste0(p_analysis, "ldsr_command.sh")

directory = p_out
if (!dir.exists(directory)) {
	dir.create(directory)
}

# CREATE MALE FEMALE COMPARISON -------------------------------------------

write("", script)

for (infile in fl_olfaction) {
	f_olfaction = paste0(p_analysis, infile)
	
	files = paste0(c(f_olfaction, files_ukbb), collapse = ",")
	# files = paste0(f_olfaction, ",", files) # add if correlation with itself is wanted
	
	p = str_match(infile, "GWASMA_(.*)_\\d+")[,2]
	f_out = str_replace(f_out_template, "PHENO", p)
	
	ldsr_command = readLines(template)
	ldsr_command = str_replace_all(ldsr_command, "FILELIST", files)
	ldsr_command = str_replace_all(ldsr_command, "OUT", paste0(p_out, f_out))
	
	write(ldsr_command, script, append = T)
	write("\n", script, append = T)
}

system(paste0("chmod +x ", script))

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
