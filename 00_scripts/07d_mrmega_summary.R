# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2025-04-02
#
# Script Description: Lookup of index variants from fixed-effect model in the MR-mega results
#
#
# Notes:
# Lookup is only performed for the best-associated phenotype
#
# pipeline_name: R10b_mrmega_summary.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
library(WriteXLS)
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

p_mrmega_results = "13_heterogeneity/mrmega/output/"
loci = fread(paste0(path_locus_definition, "locus_definition_rsID.csv"), dec = ",")
fl_mrmega = list.files(p_mrmega_results, pattern = ".gz$", full.names = TRUE)

gw_sig_threshold = 5e-8
gw_bonf_threshold = gw_sig_threshold / 13

max.cores = 30
# LOAD MRMEGA DATA --------------------------------------------------------

n.cores = min(max.cores, length(fl_mrmega))
my.cluster = parallel::makeCluster(n.cores, type = "FORK")
tmp_out = clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
d_mega = foreach(infile = fl_mrmega, .packages = "data.table") %dopar% {
  d_tmp = fread(infile, nThread = 1)
  d_tmp
}
parallel::stopCluster(cl = my.cluster)

names(d_mega) = str_match(fl_mrmega, "/mrmega_(.*)\\.result\\.gz$")[, 2]

# LOOKUP ------------------------------------------------------------------

res = loci[, .(region, markerID, phenotype, betaFEM, seFEM, pFEM)]
res = foreach(r = 1:nrow(res)) %do% {
  row = res[r]
  snp = d_mega[[row$phenotype]][MarkerName == row$markerID]
  stopifnot(nrow(snp) == 1)
  row = cbind(row, snp[, !c("MarkerName")])
}
res = rbindlist(res)
res[, gw_sig := P.value_association < 5e-8]
res[, sig_het := P.value_ancestry_het < 0.05]
res[pFEM < gw_sig_threshold, genome_wide_significance := "trait-wise"]
res[pFEM < gw_bonf_threshold, genome_wide_significance := "study-wide"]

# formatting
res = res[, c(
  "region",
  "markerID",
  "Chromosome",
  "Position",
  "phenotype",
  "betaFEM",
  "seFEM",
  "pFEM",
  "genome_wide_significance",
  "EA",
  "NEA",
  "EAF",
  "Nsample",
  "Ncohort",
  "Effects",
  "ndf_association",
  "P.value_association",
  "ndf_ancestry_het",
  "P.value_ancestry_het",
  "ndf_residual_het",
  "P.value_residual_het",
  "sig_het"
)]

WriteXLS(
  res,
  "13_heterogeneity/mrmega/index_variants_mrmega_supplement.xlsx",
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
