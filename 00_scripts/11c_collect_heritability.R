# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-03
#
# Script Description: Collect heritability information from the log-files and format them into a table
#
#
# Notes:
#
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

p_data = "19_ldsr_single_studies/heritability/"
f_out = "19_ldsr_single_studies/heritablility_overview.xlsx"

fl = list.files(p_data, pattern = "\\.log$")


# FUNCTIONS ---------------------------------------------------------------

get_value = function(d, pattern = "", parentheses = FALSE) {
  lines = d[str_detect(d, pattern)]
  stopifnot(length(lines) == 1)
  num = str_match(lines, ": ([-|\\d|.]+)")[, 2] %>% as.numeric()
  if (parentheses) {
    num = str_match(lines, "\\(([-|.|\\d]+)\\)")[, 2] %>% as.numeric()
  }

  return(num)
}

# COLLECT INFORMATION -----------------------------------------------------

results = foreach(infile = fl) %do% {
  d = readLines(paste0(p_data, infile))
  phenotype = str_match(infile, "(.+)_h2\\.log")[2]

  h2 = get_value(d, "h2:")
  h2_se = get_value(d, "h2:", parentheses = TRUE)
  lambda = get_value(d, "Lambda GC:")
  mean_chisq = get_value(d, "Mean Chi\\^2:")
  intercept = get_value(d, "Intercept:")
  intercept_se = get_value(d, "Intercept:", parentheses = TRUE)

  res = data.table(phenotype, h2, h2_se, intercept, intercept_se, lambda, mean_chisq)
}
results = rbindlist(results)

# OUTPUT ------------------------------------------------------------------

setorder(results, -intercept)
WriteXLS::WriteXLS(
  results,
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
