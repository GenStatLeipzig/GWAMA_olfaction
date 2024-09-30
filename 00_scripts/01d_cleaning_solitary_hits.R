# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2023-12-18
#
# Script Description: Cleansing the locus definition of solitary variants without support. 
# Top-SNPs that have no neighbouring SNP with at least suggestive significance are filtered
# out and a new ID is introduced to keep continous enumeration.
#
#
# Notes:
# pipeline_name: 05e_cleaning_locus_definition.R
#
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

locus_definition_file = paste0(path_locus_definition, "locus_definition_suggestive.csv") #TODO check locus definition file
file_out = str_replace(locus_definition_file, ".csv", "_cleaned.csv")

suggestive_sig = 10^-6 #TODO set threshold

max.cores = 40

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

# LOAD DATA ---------------------------------------------------------------

locus_def = fread(locus_definition_file, dec = ",")

#load the phenotype data
phenotypes_to_load = unique(locus_def$phenotype)
n.cores = min(max.cores, length(phenotypes_to_load))

my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out <- clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
pheno_data = foreach (pheno = phenotypes_to_load, .packages = "data.table") %dopar% {
  f = list.files(path_data, pattern = pheno, full.names = T)
  data = fread(f, nThread = 1)
  data = data[nWeightedMAF > maf_filter &
                nWeightedInfoScore > info_filter &
                I2 < i2_filter  &
                numberStudies >= n_studies_filter &
                numberLargeStudies >= n_large_studies_filter,]
  data
}
parallel::stopCluster(cl = my.cluster)

names(pheno_data) = phenotypes_to_load


# DEFINE REGIONS WITHOUT SUPPORT FOR TOP SNP ------------------------------

regions_to_remove = c()

for(row in 1:nrow(locus_def)){
  region = locus_def[row, region]
  reg_chrom = locus_def[row, chrom]
  reg_start = locus_def[row, region_start]
  reg_end = locus_def[row, region_end]
  reg_pheno = locus_def[row, phenotype]
  reg_snp = locus_def[row, markerID]
  
  
  data = pheno_data[[reg_pheno]]
  data = data[chrom == reg_chrom & pos >= reg_start & pos <= reg_end, ]
  data_suggestive_neigbours = data[markerID != reg_snp & pFEM <= suggestive_sig,]
  n_sugg_neigbours = nrow(data_suggestive_neigbours)
  
  if(n_sugg_neigbours == 0){
    regions_to_remove = c(regions_to_remove, region)
  }
  
}


# PROCESSING --------------------------------------------------------------

data = fread(locus_definition_file, dec = ",")

data = data[!(region %in% regions_to_remove),]

#copy old ids and create new continous ones
data[,regionID_unfiltered := region]
data[,region:= c(1:nrow(data))]

#output
write.table(
  data,
  file = file_out,
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
