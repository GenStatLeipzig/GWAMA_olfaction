# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-05-24
#
# Script Description: For regions 7,8,9 condition on the index variants of the other two loci
# to investigate whether they are independent despite their relative proximity
#
#
# Notes:
#
# pipeline_name: T8_conditioning_chr11.R

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

f_locus_def = paste0(path_locus_definition, "locus_definition.csv")
region_filter = c(7,8,9)

path_cojo_input = "07d_chr11_cojo/input/"
path_cojo_output = "07d_chr11_cojo/output/"

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

locus_definition_file = "locus_definition.csv" #TODO check file
pheno = "pineapple_all"
# COJO DATA PREPARATION ---------------------------------------------------

loci = fread(paste0(path_locus_definition, locus_definition_file), dec = ",")
loci = loci[region %in% region_filter]
region_chrom = unique(loci$chrom)

f = list.files(path_data, pheno, full.names = T)

data_tmp = fread(f, nThread = 30)

# qc filtering
data_tmp = data_tmp[nWeightedMAF > maf_filter &
                      nWeightedInfoScore > info_filter &
                      I2 < i2_filter  &
                      numberStudies >= n_studies_filter &
                      numberLargeStudies >= n_large_studies_filter,]

#note that A1 should be the effect allele and A2 the other allele
cojo_input = data_tmp[(chrom == region_chrom) & (pos >= min(loci$region_start)) & (pos <= max(loci$region_end)),
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
id_conv_file = file=paste0(path_cojo_input, "cojo_id_conversion.csv")

fwrite(
  id_conversion,
  id_conv_file,
  row.names = FALSE,
  quote = F,
  sep = ";",
  dec = ","
)
cojo_input[mismatch, SNP := new_id]    
  


# PERFORMING COJO ---------------------------------------------------------

#select sex wise index variants
snp1 = loci[1,markerID]
snp2 = loci[2,markerID]
snp3 = loci[3,markerID]

if(snp1 %in% id_conversion$id_meta){
  snp1 = id_conversion[id_meta == snp1, ld_ref_id]
}
if(snp2 %in% id_conversion$id_meta){
  snp2 = id_conversion[id_meta == snp2, ld_ref_id]
}
if(snp3 %in% id_conversion$id_meta){
  snp3 = id_conversion[id_meta == snp3, ld_ref_id]
}

snp1 = paste0("chr", snp1)
snp2 = paste0("chr", snp2)
snp3 = paste0("chr", snp3)

snps = c(snp1, snp2, snp3)
names(snps) = as.character(region_filter)

cojo_file = paste0(path_cojo_input,
                   "cojo_regions_7_8_9.ma")
fwrite(
  cojo_input,
  file = cojo_file ,
  row.names = F,
  quote = F,
  sep = " "
)

snp_file_1 = str_glue("{path_cojo_input}region{loci[1,region]}.txt")
snp_file_2 = str_glue("{path_cojo_input}region{loci[2,region]}.txt")
snp_file_3 = str_glue("{path_cojo_input}region{loci[3,region]}.txt")

writeLines(c(snp2,snp3), snp_file_1)
writeLines(c(snp1,snp3), snp_file_2)
writeLines(c(snp1,snp2), snp_file_3)

for(snp_file in c(snp_file_1, snp_file_2, snp_file_3)){
  out = str_remove(snp_file, ".txt")
  out = str_replace(out, "input", "output")
  
  gctaCommand = str_glue("{gctaCall} --bfile {path_plink}chr_{region_chrom} --chr {region_chrom} --cojo-cond {snp_file} --maf 0.01 --cojo-p 5e-8 --cojo-file {cojo_file} --out {out}")
system(gctaCommand)
  
}

result = foreach(r = region_filter) %do% {
  data_cond = fread(str_glue("{path_cojo_output}region{r}.cma.cojo"))
  snp = snps[as.character(r)]
  data_cond = data_cond[SNP == snp]
  
  cond_snps = readLines(str_glue("{path_cojo_input}region{r}.txt"))
  cond_snps = str_remove_all(paste(cond_snps, collapse = ","), "chr")
  result = data.table(region = r,
                      snp = loci[region == r, markerID],
                      cond_on = cond_snps,
                      p_ori = loci[region == r, pFEM],
                      p_cond = data_cond$pC
                      )
}

result = rbindlist(result)
result[,lost_gw_sig := p_cond>=5E-8]
fwrite(
  result,
  paste0(path_cojo_output,"conditioned_p_values.csv"),
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
