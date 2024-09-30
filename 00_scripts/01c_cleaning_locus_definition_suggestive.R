# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2023-12-18
#
# Script Description: Cleansing the locus definition with suggestive hits from 
# the genomewide sig. hits so that they can be analysed idependently analysed
#
#
# Notes: overwrites original input file
# pipeline_name: 05e_cleaning_locus_definition_suggestive.R
#

# INIT --------------------------------------------------------------------

rm(list = ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

locus_definition_file = paste0(path_locus_definition, "locus_definition_suggestive.csv") #TODO check locus def file
# file_out = str_replace(locus_definition_file, ".csv", "_cleaned.csv")
file_out = locus_definition_file

# LOAD DATA ---------------------------------------------------------------

locus_def = fread(locus_definition_file, dec = ",")


# SELECT GWSIG LOCI TO REMOVE ---------------------------------------------

regions_to_remove = locus_def[locus_def$genomewide_sig,]$region

# PROCESSING --------------------------------------------------------------

data = fread(locus_definition_file, dec = ",")

data = data[!(region %in% regions_to_remove),]

#copy old ids and create new continous ones
data[,regionID_unfiltered := region]
data[,region:= c(1:nrow(data))]
data[,region:= paste0("S", region)] #add s for suggestive

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
