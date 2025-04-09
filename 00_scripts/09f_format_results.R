# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-09
#
# Script Description: Convert .log file into csv format
#
#
# Notes:
# pipelin_name: R12c_format_results.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
library("openxlsx")

# VARIABLES ---------------------------------------------------------------

start_row = "p1"
end_row = "Analysis finished"

p_analysis = "18_ldsr_coffee/"
p_ukbb = "sum_stats_formatted/"
p_out = paste0(p_analysis, "genetic_correlation/")
f_out = paste0(p_out, "genetic_correlation.xlsx")

# PROCESSING --------------------------------------------------------------

fl = list.files(p_out, pattern = "\\.log$", full.names = TRUE)

d = foreach(f_in = fl) %do% {
	d = readLines(f_in)
	start = which(str_detect(d, start_row))
	end = which(str_detect(d, end_row)) - 2 # -2 due to the empty line above
	d = d[start:end]
	
	writeLines(d, "tmp_table.csv")
	d = fread("tmp_table.csv")
	file.remove("tmp_table.csv")
	d
}
d = rbindlist(d)

# FORMATTING --------------------------------------------------------------

d[,p1 := str_match(p1, "GWASMA_(.*)_\\d+")[,2]]
d[str_detect(p2, "_CI_"), p2 := "coffee intake"]
d[str_detect(p2, "_TI_"), p2 := "tea intake"]

WriteXLS::WriteXLS(
	d,
	f_out,
	AdjWidth = TRUE,
	AutoFilter = TRUE,
	BoldHeaderRow = TRUE
)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
