# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-11-29
#
# Script Description: Perform a sensitivity analysis to check if heterogeneity is driven by a single study.
# For each index variant meta analysis is repeated by sequentially leaving out one of the studies.
#
#
# Notes:
# pipeline_name: R01b_I2_sensitivity.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
library(openxlsx)
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

loci = fread(paste0(path_locus_definition, "locus_definition_rsID.csv"), dec = ",")

check_p_val = FALSE # check if newly calculated p-values are identical to original ones

studies = c("LIFE", "CHRIS", "Rhineland", "ARIC")

# FUNCTIONS ---------------------------------------------------------------

perform_sensitivity = function(i = 1) {
  beta_sens = beta_all[-i]
  se_sens = se_all[-i]
  study = studies[i]
  sens = purrr::map2(beta_sens, se_sens, ~ metafor::rma(yi = .x, sei = .y, method = "FE"))

  col1 = paste0("I2_wo_", study)
  col2 = paste0("het_p_wo_", study)

  res = data.table(col1 = purrr::map(sens, ~ round(.x$I2, 1)), col2 = purrr::map(sens, ~ round(.x$QEp, 4)))
  names(res) = c(col1, col2)
  return(res)
}

compare_change = function(x, y, vals = c(0, 1, -1)) {
  stopifnot(length(x) == length(y))
  result = ifelse(y > x, vals[2], vals[3])
  equal_pos = which(y == x)
  result[equal_pos] = vals[1]

  return(result)
}

# PERFORME FIXED EFFECT META ANALYSIS -------------------------------------
het_loci = loci[-7] # exclude loci that are not independent

beta_all = het_loci[, .(beta.LIFE_EUR, beta.CHRIS_EUR, beta.Rhineland_EUR, beta.ARIC_EUR)]
beta_all = transpose(beta_all)
names(beta_all) = het_loci$rsID

se_all = het_loci[, .(se.LIFE_EUR, se.CHRIS_EUR, se.Rhineland_EUR, se.ARIC_EUR)]
se_all = transpose(se_all)
names(se_all) = het_loci$rsID

all_studies = purrr::map2(beta_all, se_all, ~ metafor::rma(yi = .x, sei = .y, method = "FE"))

# check that values are identical to original meta analysis to verify calculation
results = data.table(region = het_loci$region, rsID = het_loci$rsID, phenotype = het_loci$phenotype, I2_ori = purrr::map(all_studies, ~ round(.x$I2, 1)), het_p_ori = purrr::map(all_studies, ~ round(.x$QEp, 4)))

stopifnot(all(results$I2 == het_loci$I2))

if (check_p_val) {
  hetp_values = foreach(r = 1:nrow(het_loci)) %do% {
    d_row = het_loci[r]
    infile = list.files("02_MetaGWAS/aa_excluded/", pattern = d_row$phenotype, full.names = TRUE)
    d_snp = fread(infile, nThread = 30)
    het_p = d_snp[markerID == d_row$markerID, HetPVal]
  }

  t1 = hetp_values
  t1 = unlist(t1)
  t1 = round(t1, 2)
  t2 = unlist(results$het_p)
  t2 = round(t2, 2)
  stopifnot(all.equal(t1, t2))
}

## perform sensitivity analysis ####
sens_results = foreach(i = 1:nrow(beta_all)) %do% {
  perform_sensitivity(i)
}

results = Reduce(cbind, sens_results, init = results)

results = results[, .(region, rsID, phenotype, I2_ori, I2_wo_LIFE, I2_wo_CHRIS, I2_wo_Rhineland, I2_wo_ARIC, het_p_ori, het_p_wo_LIFE, het_p_wo_CHRIS, het_p_wo_Rhineland, het_p_wo_ARIC)]

r_excel = results[, .(rsID, phenotype, I2_ori, I2_wo_LIFE, I2_wo_CHRIS, I2_wo_Rhineland, I2_wo_ARIC, het_p_ori, het_p_wo_LIFE, het_p_wo_CHRIS, het_p_wo_Rhineland, het_p_wo_ARIC)]
wb = createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", r_excel)

# Add conditional formatting
for (col in 4:7) {
  conditionalFormatting(
    wb,
    sheet = "Sheet1",
    cols = col,
    rows = 2:(nrow(results) + 1),
    rule = "<C2",
    style = createStyle(fontColour = "black", bgFill = "green4") # Green
  )

  conditionalFormatting(
    wb,
    sheet = "Sheet1",
    cols = col,
    rows = 2:(nrow(results) + 1),
    rule = ">C2",
    style = createStyle(fontColour = "black", bgFill = "darkred") # Green
  )
}
conditionalFormatting(
  wb,
  sheet = "Sheet1",
  cols = 8:12,
  rows = 2:(nrow(results) + 1),
  rule = "< 5e-2",
  style = createStyle(fontColour = "black", bgFill = "darkred") # Green
)



# Save the workbook
saveWorkbook(wb, "13_heterogeneity/sensitivity_analysis_highlight.xlsx", overwrite = TRUE)



WriteXLS::WriteXLS(
  results,
  "13_heterogeneity/sensitivity_analysis.xlsx",
  AdjWidth = TRUE,
  AutoFilter = TRUE,
  BoldHeaderRow = TRUE
)


# VISUALIZATION -----------------------------------------------------------
i2_cols = names(results)[str_detect(names(results), "I2")]
plot_data = results[, ..i2_cols]
plot_data[, (i2_cols) := lapply(.SD, unlist), .SDcols = i2_cols]
# cols = as.matrix(plot_data)
change = purrr::modify(plot_data, ~ compare_change(plot_data$I2_ori, .x))

rnames = paste("Locus", results$region, "-", results$rsID)
cnames = str_replace_all(names(change), "_", " ") %>%
  str_replace_all(., "ori", "meta") %>%
  str_replace_all(., "wo ", "-")

pdf("13_heterogeneity/sensitivity.pdf")
pheatmap(change,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  display_numbers = plot_data,
  color = c("#BF616A", "white", "#A3BE8C"),
  legend = FALSE,
  labels_row = rnames,
  labels_col = cnames,
  gaps_col = 1
)
dev.off()



# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
