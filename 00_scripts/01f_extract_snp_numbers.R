# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-07-29
#
# Script Description: Correct locus definition case numbers
#
#
# Notes: The locus definition script does not correctly determine the number of 
# significant SNPs at certain loci (when no merge process is needed). NR_SNPs only
# counts the number of filtered SNPs (0 in case no merge process happens)
# This script replaces the NR_SNPs value with the correct count of QC SNPs 
# below the significance threshold within the region boundries
#
# pipeline_name: 05f_extract_snp_numbers.R

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------
f_loci = paste0(path_locus_definition, "locus_definition.csv") #TODO check locus definition file
f_out = str_replace(f_loci, "\\.csv$", "\\_snp_numbers.csv")

loci = fread(f_loci, dec = ",")

#TODO check significance threshold
p_filter = 5e-8
# p_filter = 1e-6 #use for suggestive locus definition

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

max.cores = 40

# LOAD PHENOTYPES ---------------------------------------------------------

f_phenos = list.files(path_data)

n.cores = min(max.cores, length(f_phenos))
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out <- clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
data = foreach(f = f_phenos, .packages = "data.table") %dopar%{
  dat = fread(paste0(path_data, f), nThread = 1)
  
  #quality filter
  dat = dat[nWeightedMAF > maf_filter &
              nWeightedInfoScore > info_filter &
              I2 < i2_filter  &
              numberStudies >= n_studies_filter &
              numberLargeStudies >= n_large_studies_filter,]
}
parallel::stopCluster(cl = my.cluster)

names(data) = str_match(f_phenos, "GWASMA_([^_]+_[^_]+)_")[,2]


# EXTRACT LOCI NUMBERS ----------------------------------------------------
# ori_loci = copy(loci)
for(r in 1:nrow(loci)){
  d = data[[loci[r, phenotype]]]
  d = d[chrom == loci[r, chrom] &
          pos >= loci[r, region_start] &
          pos <= loci[r, region_end] & 
          pFEM < p_filter]
  loci[r, NR_SNPs := nrow(d)]
}

write.table(
  loci,
  file = f_out,
  col.names = T,
  row.names = F,
  quote = F,
  sep = ";",
  dec = ","
)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
