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
# pipeline_name: R06c_format_results.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
library(openxlsx)

# VARIABLES ---------------------------------------------------------------

start_row = "p1"
end_row = "Analysis finished"

p_analysis = "15_ldsr_diseases/"
p_ukbb = "UKBB_phenotypes/"
p_out = paste0(p_analysis, "genetic_correlation/")
f_in = paste0(p_out, "rg_score_all_diseases.log")
f_phenoslct = paste0("additional_information/", "ldsc_selected_phenotypes.xlsx")
f_out = str_replace(f_in, "\\.log$", "\\.xlsx")



# PROCESSING --------------------------------------------------------------

d = readLines(f_in)
start = which(str_detect(d, start_row))
end = which(str_detect(d, end_row)) - 2 # -2 due to the empty line above
d = d[start:end]

writeLines(d, "tmp_table.csv")
d = fread("tmp_table.csv")
file.remove("tmp_table.csv")


# FORMATTING --------------------------------------------------------------

d[, p1 := "score_all"]
d[, p2 := str_split_i(p2, "__", 1)]
d[, p2 := str_split_i(p2, "//", 2)]
d[1, p2 := "score_all"]

WriteXLS::WriteXLS(
  d,
  paste0(f_out),
  AdjWidth = TRUE,
  AutoFilter = TRUE,
  BoldHeaderRow = TRUE
)


# COMBINE WITH STUDY INFORMATION ------------------------------------------


d_ukbb = read.xlsx(f_phenoslct) %>% as.data.table()

# format descriptors to match on LDSC results
matching_description = ifelse(
  str_starts(d_ukbb$description, "Non-cancer illness code"),
  str_split_i(d_ukbb$description, ": ", 2),
  d_ukbb$description
) %>%
  tolower() %>%
  str_replace_all(., " ", "_")
d_ukbb[, matching_description := matching_description]
stopifnot(all(d_ukbb$matching_description %in% d$p2))

d_ukbb = d_ukbb[, .(phenotype, description, source, sex, ldsc_sumstat_dropbox, matching_description)]

d = merge(d, d_ukbb, by.x = "p2", by.y = "matching_description", all.x = TRUE, sort = FALSE)
cols = names(d)
cols[1] = "p1"
cols[2] = "p2"
d = d[, ..cols]
setorder(d, -rg, na.last = TRUE)

WriteXLS::WriteXLS(
  d,
  str_replace(f_out, "\\.xlsx$", "_supplement.xlsx"),
  AdjWidth = TRUE,
  AutoFilter = TRUE,
  BoldHeaderRow = TRUE,
  AllText = TRUE
)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
