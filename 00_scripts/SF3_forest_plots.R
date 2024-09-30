# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2023-09-27
#
# Script Description: Creates forest plots for the significant loci in the locus definition
#
#
# Notes:
# run only on forostar, throws issues on aman due to different library versions
# RsID is used as an identifier
# individual plots can be merged with PdfMerge
# pipeline_name: 07b_forest_plots_publication.R

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server

suppressPackageStartupMessages(library(forestplot))
suppressPackageStartupMessages(library(dplyr))

setwd(projectpath)


# VARIABLES ---------------------------------------------------------------
n.cores = 40 #TODO check number of needed cores

locus_definition_file = "locus_definition_rsID.csv" #TODO check locus definition
region_data = fread(paste0(path_locus_definition, locus_definition_file), dec = ",")

# DEFINE PROCEDURE --------------------------------------------------------

plot_forest = function(row){
  region_id = region_data[row, region]
  snp = region_data[row, markerID]
  rsID = region_data[row, rsID]
  phenotypes = str_split(region_data[row, phenotypes_in_region], " \\| ")[[1]]
  
  for (phenotype in phenotypes){
    meta_analysis_file = list.files(path_data, pattern = paste0("_", phenotype, ".*\\.gz$"))
    data_meta = fread(paste0(path_data, meta_analysis_file), nThread = 1)
    
    #TODO add/remove ARIC AA
    row_data = data_meta[markerID == snp,]
    
    plot_data = data.frame(
      mean = c(row_data$beta.LIFE_EUR, 
               row_data$beta.Rhineland_EUR,
               row_data$beta.ARIC_EUR,
               # row_data$beta.ARIC_AFR,
               row_data$beta.CHRIS_EUR,
               row_data$betaFEM),
      
      lower = c(row_data$beta.LIFE_EUR-1.96*row_data$se.LIFE_EUR,
                row_data$beta.Rhineland_EUR-1.96*row_data$se.Rhineland_EUR,
                row_data$beta.ARIC_EUR-1.96*row_data$se.ARIC_EUR,
                # row_data$beta.ARIC_AFR-1.96*row_data$se.ARIC_AFR,
                row_data$beta.CHRIS_EUR-1.96*row_data$se.CHRIS_EUR,
                row_data$betaFEM-1.96*row_data$seFEM),
      
      upper = c(row_data$beta.LIFE_EUR+1.96*row_data$se.LIFE_EUR,
                row_data$beta.Rhineland_EUR+1.96*row_data$se.Rhineland_EUR,
                row_data$beta.ARIC_EUR+1.96*row_data$se.ARIC_EUR,
                # row_data$beta.ARIC_AFR+1.96*row_data$se.ARIC_AFR,
                row_data$beta.CHRIS_EUR+1.96*row_data$se.CHRIS_EUR,
                row_data$betaFEM+1.96*row_data$seFEM),
      
      study = c("LIFE-Adult",
                "Rhineland",
                "ARIC EA",
                # "ARIC AA",
                "CHRIS",
                "Meta-Result"),
      
      N = c(row_data$n.LIFE_EUR,
            row_data$n.Rhineland_EUR,
            row_data$n.ARIC_EUR,
            # row_data$n.ARIC_AFR,
            row_data$n.CHRIS_EUR,
            row_data$totalN),
      
      is.summary = c(F,
                     F,
                     F,
                     # F,
                     F,
                     T)
    )
    
    header = data.frame(mean = NA, lower = NA, upper = NA, study = "Study", N = "N", is.summary = T)
    plot_data = rbind(header, plot_data)
    
    pdf(paste0("04_mh_qq_plots_overview_stats/forest_plots_region_", region_id,"_",phenotype, ".pdf"), width = 7, height = 5)
    plot = forestplot(plot_data,
               labeltext = c(study, N),
               title = str_glue("{rsID}\nLocus: {region_id}, phenotype: {str_replace(str_to_lower(phenotype), '_', ' ')}, I2: {round(row_data$I2,3)}"),
               xlab = "Effect",
               is.summary = is.summary,
               new_page = F,
               # lwd.ci = 1,
               lwd.zero = 3,
               col = fpColors(lines = "black")
               )
    print(plot)
    dev.off()
  }
 
}


# PARALLLEL PROCESSING ----------------------------------------------------

my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out <-clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
result = foreach(
  row = c(1:nrow(region_data)),
  .packages = c("data.table", "stringr", "forestplot", "dplyr")
) %dopar% {
  plot_forest(row)
}
parallel::stopCluster(cl = my.cluster)


# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
