# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-02-09
#
# Script Description: Performs COJO conditional on known literature hits
# to identify possible secondary hits
#
#
# Notes:
# pipeline_name: 14c_identify_secondary_hits.R
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server


# VARIABLES ---------------------------------------------------------------

novelty_file = "additional_information/novelty_status.csv" #manual table of intersection between defined loci and 
locus_def_file = "locus_definition.csv" #TODO check locus definition
path_cojo_input = "07b_cojo_cond_lit/input/"
path_cojo_output = "07b_cojo_cond_lit/output/"


# LOAD DATA ---------------------------------------------------------------

novelty = fread(novelty_file)
novelty = novelty[!is.na(lies_within_regions)]
novelty[,chrom := as.integer(chrom)]

locus_def = fread(paste0(path_locus_definition, locus_def_file), dec = ",")

#QC
p_filter = 5 * 10 ^ -8
maf_filter = 0.01
info_filter = 0.8 #TODO default was 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

# PROCESSING --------------------------------------------------------------

for (row in 1:nrow(novelty)){
  
  lit_chrom = novelty[row, chrom]
  lit_pos = novelty[row, pos]
  lit_region = novelty[row, lies_within_regions]
  
  # load data for replication
  region_data = locus_def[region==lit_region,]
  phenotypes = str_split(region_data$phenotypes_in_region, " \\| ", simplify = T)
  region_start = region_data$region_start
  region_end = region_data$region_end
  
  # for each phenotype load the snp data, perform qc-filtering and condition on lit snp
  
  for (phenotype in phenotypes){
    pheno_file = list.files(path_data, pattern = phenotype, full.names = T)
    data_tmp = fread(pheno_file, nThread = 20)
    
    #selcect lit snp to conditon on 
    lit_snp = data_tmp[chrom==lit_chrom & pos==lit_pos]$markerID
    lit_snp = paste0("chr", lit_snp)
    if(length(lit_snp)!=1){
      message("Warning: SNP not found or not unique for: ", novelty[row,rsID], "\n")
    }
    
    #check that snp is compatible with LD reference (some IDs have changed)
    bim = fread(str_glue("{path_plink}chr_{lit_chrom}.bim"), nThread = 10)
    if(!(lit_snp %in% bim$V2)){
      id_parts = str_split(lit_snp, ":", simplify = T)
      lit_snp = paste(id_parts[1], id_parts[2], id_parts[4], id_parts[3], sep=":")
      data_tmp[chrom==lit_chrom & pos==lit_pos, markerID := str_remove(lit_snp, "chr")]
      if(!(lit_snp %in% bim$V2)){
        next
      }
    }
    
    #qc and select regional snps
    # no qc because lit snp might be filtered out otherwise
    # data_tmp = data_tmp[nWeightedMAF > maf_filter &
    #                       nWeightedInfoScore > info_filter &
    #                       I2 < i2_filter  &
    #                       numberStudies >= n_studies_filter &
    #                       numberLargeStudies >= n_large_studies_filter,]
    
    #note that A1 should be the effect allele and A2 the other allele
    cojo_input = data_tmp[(chrom == lit_chrom) & (pos >= region_start) & (pos <= region_end),
                          c("markerID", "ea", "aa", "nWeightedEAF", "betaFEM", "seFEM", "pFEM", "totalN")]
    # modify ID to work with plink reference files
    cojo_input[,markerID := paste0("chr", markerID)]
    
    setnames(cojo_input, c("SNP", "A1", "A2", "freq", "b", "se", "p", "N"))
    
    # some IDs are switched between the LIFE reference and COJO -> change the IDs
    # for the mismatching ones
    bim = fread(str_glue("{path_plink}chr_{lit_chrom}.bim"), nThread = 10)
    
    mismatch = !(cojo_input$SNP %in% bim$V2) #mismatches with original ID
    
    mismatching_snps = cojo_input[mismatch]$SNP
    new_id = transpose(str_split(mismatching_snps, ":"))
    new_id = paste(new_id[[1]], new_id[[2]], new_id[[4]], new_id[[3]], sep = ":")
    
    id_conversion = data.table(id_meta = cojo_input[mismatch]$SNP, ld_ref_id = new_id)
    id_conv_file = paste0(path_cojo_input, "cojo_id_conversion_", lit_region, "_", phenotype, "__chr", lit_chrom,".csv")
    
    fwrite(
      id_conversion,
      id_conv_file,
      row.names = FALSE,
      quote = F,
      sep = ";",
      dec = ","
    )
    cojo_input[mismatch, SNP := new_id]    
    
    
    # output
    snp_file = paste0(path_cojo_input, "cojo_input_region_", lit_region, "_", phenotype, "__chr", lit_chrom,".ma")
    fwrite(cojo_input, file=snp_file, row.names = F, quote = F, sep = " ")
    
    #write snp file with lit snp to condition on 
    lit_snp_file = paste0(path_cojo_input, "cojo_lit_snp_", lit_region, "_", phenotype, "__chr", lit_chrom,".ma")
    
    fileConn = file(lit_snp_file)
    writeLines(lit_snp, fileConn)
    close(fileConn)
    
    #perform cojo cond
    output_file = str_replace_all(snp_file, "input", "output")
    output_file = str_remove(output_file, "\\.ma")
    
    
    gctaCommand = str_glue("{gctaCall} --bfile {path_plink}chr_{lit_chrom} --chr {lit_chrom} --cojo-cond {lit_snp_file} --maf 0.01 --cojo-p 5e-8 --cojo-file {snp_file} --out {output_file}")
    system(gctaCommand)
    
    
    #collect additional gw hits
    cond_data = fread(paste0(output_file, ".cma.cojo"))
    cond_data = cond_data[pC < p_filter]
    
    fwrite(
      cond_data,
      file = paste0(path_cojo_output, "secondary_hits_", lit_region, "_", phenotype, ".csv"),
      row.names = FALSE,
      quote = F,
      sep = ";",
      dec = ","
    )
  }
}

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
