# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2025-01-17
#
# Script Description: For the colocalization analyses, genetic data of the regions is needed
# As the wholde data set is to large, we crop out the regions of interest.
#
#
# Notes:
# pipeline_name: tmp_david_coloc_data.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

loci = fread("05_locus_definition/locus_definition_rsID.csv", dec = ",")
loci = loci[!region == 7]

phenos = str_split(loci$phenotypes_in_region, " \\| ")
names(phenos) = loci$region

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

# LOAD DATA ---------------------------------------------------------------

phenos_to_load = unlist(phenos) %>% unique()

n.cores = min(30, length(phenos_to_load))
my.cluster = parallel::makeCluster(n.cores, type = "FORK")
tmp_out = clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
d = foreach(p = phenos_to_load, .packages = loadedNamespaces()) %dopar% {
  infile = list.files(path_data, pattern = p, full.names = TRUE)
  stopifnot(length(infile) == 1)
  d_tmp = fread(infile, nThread = 1)
  d_tmp = d_tmp[nWeightedMAF > maf_filter &
    nWeightedInfoScore > info_filter &
    I2 < i2_filter &
    numberStudies >= n_studies_filter &
    numberLargeStudies >= n_large_studies_filter, ]
}
parallel::stopCluster(cl = my.cluster)
names(d) = phenos_to_load


# PREPARE REGION DATA -----------------------------------------------------

directory = "tmp_data_for_coloc/"
if (!dir.exists(directory)) {
  dir.create(directory)
}


d_coloc = foreach(r = 1:nrow(loci)) %do% {
  row = loci[r]
  d_region = foreach(p = phenos[[as.character(row$region)]]) %do% {
    d_pheno = d[[p]]
    d_pheno = d_pheno[chrom == row$chrom & pos >= row$region_start & pos <= row$region_end]
    d_pheno[, region := row$region]
    d_pheno[, phenotype := p]
    is_top_pheno = p == row$phenotype
    d_pheno[, is_top_phenotype := is_top_pheno]
  }
  d_region = rbindlist(d_region)
}
d_coloc = rbindlist(d_coloc)


# OUTPUT ------------------------------------------------------------------

fwrite(d_coloc, "tmp_data_for_coloc/coloc_data_cropped_regions_qc_filtered.csv.gz", nThread = 30)



# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
