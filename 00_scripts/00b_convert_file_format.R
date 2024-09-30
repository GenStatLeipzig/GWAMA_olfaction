# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-09-10
#
# Script Description: Convert files from the SSF format of the GWAS catalog into 
# the internal format used for the following analyses.
#
#
# Notes:
#
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta_functional.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

max.cores = 13

# FUCTIONS ----------------------------------------------------------------

prepare_pipeline_format = function(f_gwcat){
  
  data = fread(paste0(path_raw_data, f_gwcat), na.strings = c("#NA"), nThread = 1)
  
  setnames(data, 
           c("chromosome", "base_pair_location", "effect_allele", "other_allele", 
             "beta", "standard_error", "effect_allele_frequency", "p_value", 
             "variant_id", "info", "n", "n_stud", "I2", "direction", "beta_LIFE_EUR", 
             "standard_error_LIFE_EUR", "p_value_LIFE_EUR", "n_LIFE_EUR", 
             "effect_allele_frequency_LIFE_EUR", "info_LIFE_EUR", "beta_Rhineland_EUR", 
             "standard_error_Rhineland_EUR", "p_value_Rhineland_EUR", "n_Rhineland_EUR", 
             "effect_allele_frequency_Rhineland_EUR", "info_Rhineland_EUR", 
             "beta_ARIC_EUR", "standard_error_ARIC_EUR", "p_value_ARIC_EUR", 
             "n_ARIC_EUR", "effect_allele_frequency_ARIC_EUR", "info_ARIC_EUR", 
             "beta_CHRIS_EUR", "standard_error_CHRIS_EUR", "p_value_CHRIS_EUR", 
             "n_CHRIS_EUR", "effect_allele_frequency_CHRIS_EUR", "info_CHRIS_EUR", 
             "het_chi_sq", "het_df", "het_p_value" 
           ),
           c("chrom", "pos", "ea", "aa", 
             "betaFEM", "seFEM", "nWeightedEAF", "pFEM", 
             "markerID", "nWeightedInfoScore", "totalN", "numberStudies", "I2", "direction", "beta.LIFE_EUR", 
             "se.LIFE_EUR", "p.LIFE_EUR", "n.LIFE_EUR", 
             "eaf.LIFE_EUR", "infoscore.LIFE_EUR", "beta.Rhineland_EUR", 
             "se.Rhineland_EUR", "p.Rhineland_EUR", "n.Rhineland_EUR", 
             "eaf.Rhineland_EUR", "infoscore.Rhineland_EUR", 
             "beta.ARIC_EUR", "se.ARIC_EUR", "p.ARIC_EUR", 
             "n.ARIC_EUR", "eaf.ARIC_EUR", "infoscore.ARIC_EUR", 
             "beta.CHRIS_EUR", "se.CHRIS_EUR", "p.CHRIS_EUR", 
             "n.CHRIS_EUR", "eaf.CHRIS_EUR", "infoscore.CHRIS_EUR",
             "HetChiSq", "HetDf", "HetPVal"
           ))
  
  data[, nWeightedMAF := ifelse(nWeightedEAF < 0.5, nWeightedEAF, 1-nWeightedEAF)]
  data[, markerID := str_replace_all(markerID, "_", ":")]
  data[, numberLargeStudies := numberStudies] # dummy column that needs to be present for compatibility
  data[, maf.LIFE_EUR := ifelse(eaf.LIFE_EUR < 0.5, eaf.LIFE_EUR, 1-eaf.LIFE_EUR)]
  data[, maf.Rhineland_EUR := ifelse(eaf.Rhineland_EUR < 0.5, eaf.Rhineland_EUR, 1-eaf.Rhineland_EUR)]
  data[, maf.ARIC_EUR := ifelse(eaf.ARIC_EUR < 0.5, eaf.ARIC_EUR, 1-eaf.ARIC_EUR)]
  data[, maf.CHRIS_EUR := ifelse(eaf.CHRIS_EUR < 0.5, eaf.CHRIS_EUR, 1-eaf.CHRIS_EUR)]
  
  # unused dummy cols
  
  # restore order 
  data = data[, c("markerID", "ea", "aa", "nWeightedEAF", "betaFEM", "seFEM", 
                  "pFEM", "direction", "I2", "HetChiSq", "HetDf", "HetPVal", "beta.LIFE_EUR", 
                  "se.LIFE_EUR", "p.LIFE_EUR", "n.LIFE_EUR", "eaf.LIFE_EUR", "maf.LIFE_EUR", 
                  "infoscore.LIFE_EUR", "beta.Rhineland_EUR", "se.Rhineland_EUR", 
                  "p.Rhineland_EUR", "n.Rhineland_EUR", "eaf.Rhineland_EUR", "maf.Rhineland_EUR", 
                  "infoscore.Rhineland_EUR", "beta.ARIC_EUR", "se.ARIC_EUR", "p.ARIC_EUR", 
                  "n.ARIC_EUR", "eaf.ARIC_EUR", "maf.ARIC_EUR", "infoscore.ARIC_EUR", 
                  "n.CHRIS_EUR", "maf.CHRIS_EUR", "infoscore.CHRIS_EUR", "eaf.CHRIS_EUR", 
                  "p.CHRIS_EUR", "nWeightedMAF", "numberStudies", "numberLargeStudies", 
                  "totalN", "chrom", "pos", "nWeightedInfoScore", "se.CHRIS_EUR", 
                  "beta.CHRIS_EUR")]
  
  fn_out = str_replace(f_gwcat, ".tsv.gz", paste0("_", Sys.Date(), ".gz"))
  
  fwrite(
    data,
    paste0(path_data, fn_out),
    row.names = FALSE,
    quote = FALSE,
    nThread = 1
  )

}


# PROCESSING --------------------------------------------------------------
file_list = list.files(path_raw_data, pattern = ".gz")

n.cores = min(max.cores, length(file_list))
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out <- clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
result = foreach(f_in = file_list, .packages = loadedNamespaces()) %dopar% {
  prepare_pipeline_format(f_in)
}
parallel::stopCluster(cl = my.cluster)


# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
