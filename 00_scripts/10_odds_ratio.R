# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-12
#
# Script Description: Calculate odds ratios belonging to our effects. Effect directions are harmonized beforehand to allow for easier comparison of OR sizes
#
#
# Notes:
# pipeline_name: R07_odds_ratio.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

d_loci = fread("05_locus_definition/locus_definition_rsID.csv", dec = ",")
p_out = "16_odds_ratio/"

directory = p_out
if (!dir.exists(directory)) {
  dir.create(directory)
}

gw_sig = 5e-8
gw_bonf = gw_sig/13 # Bonferroni adjusted genome-wide significance threshold

# PROCESSING --------------------------------------------------------------

d_or = d_loci[, .(region, rsID, phenotype, nWeightedMAF, pFEM, betaFEM, seFEM, aa, ea)]
d_or[, beta_harmonized := ifelse(betaFEM < 0, -betaFEM, betaFEM)]
d_or[, risk_allel := ifelse(betaFEM < 0, aa, ea)]
d_or[, or := round(exp(beta_harmonized),2)]
d_or = d_or[, .(region, rsID, phenotype, nWeightedMAF, risk_allel,pFEM, betaFEM, beta_harmonized,seFEM, or)]
names(d_or) = c("region", "rsID", "phenotype", "nWeightedMAF", "risk_allel", "pvalue", "beta_ori", "beta_harmonized", "se_beta", "OR")

d_or[, lower_ci95_or := exp(beta_harmonized - (1.96 * se_beta))]
d_or[, upper_ci95_or := exp(beta_harmonized + (1.96 * se_beta))]
d_or[, lower_ci95_or := round(lower_ci95_or, 2)]
d_or[, upper_ci95_or := round(upper_ci95_or, 2)]
d_or[pvalue < gw_sig, genome_wide_significance := "trait-wise"]
d_or[pvalue < gw_bonf, genome_wide_significance := "study-wide"]
WriteXLS::WriteXLS(
  d_or,
  paste0(p_out, "odds_ratio_publication.xlsx"),
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
