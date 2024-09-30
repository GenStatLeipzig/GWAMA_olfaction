# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-05-24
#
# Script Description: For loci with different signals in males and females perform conditioning on the sex-top SNPs to check wheter they keep their significance
#
#
# Notes:
# pipeline_name: T4_sex_conditioning.R
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------
f_locus_def = paste0(path_locus_definition, "locus_definition.csv")
l_sig_diff = c(6,11) #loci with stronges support in H3

path_cojo_input = "07c_sex_conditioning/input/"
path_cojo_output = "07c_sex_conditioning/output/"

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1


# COJO DATA PREPARATION ---------------------------------------------------

prepare_cojo_data = function(f){
  data_tmp = fread(f, nThread = 30)
  
  # qc filtering
  data_tmp = data_tmp[nWeightedMAF > maf_filter &
                        nWeightedInfoScore > info_filter &
                        I2 < i2_filter  &
                        numberStudies >= n_studies_filter &
                        numberLargeStudies >= n_large_studies_filter,]
  
  #note that A1 should be the effect allele and A2 the other allele
  cojo_input = data_tmp[(chrom == region_chrom) & (pos >= region_start) & (pos <= region_end),
                        c("markerID", "ea", "aa", "nWeightedEAF", "betaFEM", "seFEM", "pFEM", "totalN")]
  
  # modify ID to work with plink reference files
  cojo_input[,markerID := paste0("chr", markerID)]
  
  setnames(cojo_input, c("SNP", "A1", "A2", "freq", "b", "se", "p", "N"))
  
  # some IDs are switched between the LIFE reference and COJO -> change the IDs
  # for the mismatching ones
  bim = fread(str_glue("{path_plink}chr_{region_chrom}.bim"), nThread = 10)
  
  mismatch = !(cojo_input$SNP %in% bim$V2) #mismatches with original ID
  
  mismatching_snps = cojo_input[mismatch]$SNP
  new_id = transpose(str_split(mismatching_snps, ":"))
  new_id = paste(new_id[[1]], new_id[[2]], new_id[[4]], new_id[[3]], sep = ":")
  
  id_conversion = data.table(id_meta = cojo_input[mismatch]$SNP, ld_ref_id = new_id)
  pheno = str_match(f, "GWASMA_([^_]+_[^_]+)_")[2]
  id_conv_file = file=paste0(path_cojo_input, "cojo_id_conversion", region, "_", pheno, "_pos_", region_pos, "__chr", region_chrom,".csv")
  
  fwrite(
    id_conversion,
    id_conv_file,
    row.names = FALSE,
    quote = F,
    sep = ";",
    dec = ","
  )
  cojo_input[mismatch, SNP := new_id]    
  
  #output
  return(cojo_input)
  # fwrite(cojo_input, file=paste0(path_cojo_input, "cojo_input_region_", region, "_", phenotype, "_pos_", region_pos, "__chr", region_chrom,".ma"), row.names = F, quote = F, sep = " ")
}

# PERFORMING COJO ---------------------------------------------------------
regions = fread(f_locus_def, dec = ",")
regions = regions[region %in% l_sig_diff]

result = foreach(row = 1:nrow(regions)) %do% {
  phenotype = regions[row, phenotype]
  phenotype = str_split(phenotype, "_", simplify = T)[[1]]
  region_chrom = regions[row, chrom]
  region = regions[row, region]
  region_start = regions[row, region_start]
  region_end = regions[row, region_end]
  region_pos = regions[row, pos]
  
  f_m = list.files(path_data, paste0(phenotype, "_male"), full.names = T)
  f_f = list.files(path_data, paste0(phenotype, "_female"), full.names = T)
  
  cojo_in_m = prepare_cojo_data(f_m)
  cojo_in_f = prepare_cojo_data(f_f)
  
  #select sex wise index variants
  snp_m = cojo_in_m[p == min(cojo_in_m$p), SNP]
  snp_f = cojo_in_f[p == min(cojo_in_f$p), SNP]
  p_ori_m = min(cojo_in_m$p)
  p_ori_f = min(cojo_in_f$p)
  
  #output
  cojo_file_m = paste0(
    path_cojo_input,
    "cojo_region_",
    region,
    "_",
    paste0(phenotype, "_male"),
    "_pos_",
    region_pos,
    "__chr",
    region_chrom,
    ".ma"
  )
  cojo_file_f = paste0(
    path_cojo_input,
    "cojo_region_",
    region,
    "_",
    paste0(phenotype, "_female"),
    "_pos_",
    region_pos,
    "__chr",
    region_chrom,
    ".ma"
  )
  
  snp_file_m = paste0(path_cojo_input, "index_snp_region_", region, "_male.txt")
  snp_file_f = paste0(path_cojo_input, "index_snp_region_", region, "_female.txt")
  
  fwrite(
    cojo_in_m,
    file = cojo_file_m ,
    row.names = F,
    quote = F,
    sep = " "
  )
  
  fwrite(
    cojo_in_f,
    file = cojo_file_f,
    row.names = F,
    quote = F,
    sep = " "
  )
  
  
  write(snp_m, snp_file_m)
  write(snp_f, snp_file_f)
  
  
  # run the cojo command
  chrom = region_chrom
  out_m = str_remove(cojo_file_m, '\\.ma')
  out_m = str_replace(out_m, "input", "output")
  out_f = str_remove(cojo_file_f, '\\.ma')
  out_f = str_replace(out_f, "input", "output")
  
  gctaCommand = str_glue("{gctaCall} --bfile {path_plink}chr_{chrom} --chr {chrom} --cojo-cond {snp_file_f} --maf 0.01 --cojo-p 5e-8 --cojo-file {cojo_file_m} --out {out_m}")
  system(gctaCommand)
  gctaCommand = str_glue("{gctaCall} --bfile {path_plink}chr_{chrom} --chr {chrom} --cojo-cond {snp_file_m} --maf 0.01 --cojo-p 5e-8 --cojo-file {cojo_file_f} --out {out_f}")
  system(gctaCommand)
  
  
  data_m = fread(paste0(out_m, ".cma.cojo"))
  data_m = data_m[SNP == snp_m]
  data_f = fread(paste0(out_f, ".cma.cojo"))
  data_f = data_f[SNP == snp_f]
  
  result = data.table(region = rep(region, 2), sex = c("male", "female"), SNP = c(snp_m, snp_f), cond_on = c(snp_f, snp_m), p_ori = c(p_ori_m, p_ori_f), p_cond = c(data_m$pC, data_f$pC))
  
}
result = rbindlist(result)
result[,diff := p_ori - p_cond]
fwrite(
  result,
  paste0(path_cojo_output,"results_sex_conditioning.csv"),
  row.names = FALSE,
  quote = F,
  sep = ";",
  dec = ","
)


# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
