# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2023-07-21
#
# Script Description: finds the rsIDs for the loci in the locus definition
# works on matching chrom and bp position. Alleles have to be checked manually
#
# Notes: Only runs in Rstudio directly, due to the manual console input needed 
# for selecting SNPs
# pipeline_name: 05b_attatch_rs_id.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")


# variables ---------------------------------------------------------------

file_locus_definition = "locus_definition.csv" #TODO
output_file = str_replace(file_locus_definition, ".csv", "_rsID.csv")

# load_data ---------------------------------------------------------------
message("\n--------------------------\n")
message("Loading data...\n")


locus_definition = fread(paste0(path_locus_definition, file_locus_definition), dec = ",")

chroms_to_load = unique(locus_definition$chrom)

n.cores = min(40, length(chroms_to_load))
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
snp_annotation = foreach(
  chrom = chroms_to_load,
  .packages = c("data.table", "stringr"),
  .export = c("path_snp_annotation")
) %dopar% {
  if (chrom == 23)
    chrom = "X"
  annotation = fread(str_glue("{path_snp_annotation}homo_sapiens-chr{chrom}.vcf.gz"),
                     nThread = 1)
  annotation
}
parallel::stopCluster(cl = my.cluster)

names(snp_annotation) = chroms_to_load


# get rsID ----------------------------------------------------------------
rsIDs = c()
for (region in 1:nrow(locus_definition)) {
  chr = locus_definition[region, chrom]
  pos = locus_definition[region, pos]
  allele1 = str_to_upper(locus_definition[region, aa])
  allele2 = str_to_upper(locus_definition[region, ea])
  
  chrom_data = snp_annotation[[as.character(chr)]]
  snp_data = chrom_data[POS == pos, ]
  
  if (nrow(snp_data) > 1) {
    print(str_glue("Alleles: {allele1}/{allele2}"))
    print(snp_data)
    
    print("Multiple matching entries. Select row number:")
    row_number = readline()
    snp_data = snp_data[as.numeric(row_number), ]
  }
  
  if (nrow(snp_data) == 0) {
    rsID = ""
  } else {
    if ((snp_data$REF == allele1 && snp_data$ALT == allele2) ||
        (snp_data$REF == allele2 && snp_data$ALT == allele1)) {
      rsID = snp_data$ID
    } else {
      print(
        str_glue("Uncertain Alleles: We have: {allele1}/{allele2}, reference has: {snp_data$REF}/{snp_data$ALT}. Keep? (y/_)")
      )
      userInput = readline()
      if (userInput == "y") {
        rsID = snp_data$ID
      } else {
        rsID = ""
      }
    }
  }
  rsIDs = c(rsIDs, rsID)
}


# output ------------------------------------------------------------------
locus_definition[, rsID := rsIDs]

write.table(
  locus_definition,
  paste0(path_locus_definition, output_file),
  col.names = T,
  row.names = F,
  quote = F,
  sep = ";",
  dec = ","
)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished matching of rs-ID.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
