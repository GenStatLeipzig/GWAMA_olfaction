# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2023-06-13
#
# Script Description: Perform the COJO analysis for the best associated phenotype for each locus
#
#
# Notes: Performs COJO-Select analysis on the loci in the locus definition
# Script adjusts for the ID mismatches between LIFE LD reference and Meta-analysis
# pipeline_name: 09_cojo_correct_ID.R

# INIT --------------------------------------------------------------------
rm(list=ls())

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
setwd(projectpath)

locus_definition_file = "locus_definition.csv" #TODO check locus definition file
path_cojo_input = "07a_cojo_loci/input/"
path_cojo_output = "07a_cojo_loci/output/"

n.cores = 40

# qc filters
maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

# load data for COJO input creation ---------------------------------------

regions = fread(paste0(path_locus_definition, locus_definition_file), dec = ",")
# regions = regions[11,]


pheno_files_list = list.files(path_data, ".gz")
# pheno_files_list = pheno_files_list[1]

read.parallel = function(file) {
  dat = fread(paste0(path_data, file), nThread = 1)
  return(dat)
}
my.cluster <- parallel::makeCluster(n.cores,
                                    type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

#parallel processing
message("\n--------------------------\n")
message("loading SNP data...\n")
phenotype_data = foreach (file = pheno_files_list,
                          .packages = c("data.table")) %dopar% {
                            read.parallel(file)
                          }
parallel::stopCluster(cl = my.cluster)

getPheno = function(f) {
  phenotype = str_match(f, "GWASMA_(\\D+)_\\d")[2]
  return(phenotype)
}

names(phenotype_data) = sapply(pheno_files_list, getPheno)


# filter SNPs from defined regions as COJO input --------------------------
message("\n--------------------------\n")
message("Preparing COJO input files...\n")
loci_list = c()
for(row in 1:nrow(regions)){
  
  phenotypes = regions[row, phenotypes_in_region]
  phenotypes = str_split(phenotypes, " \\| ", simplify = T)
  # phenotype = regions[row, phenotype]
  region_chrom = regions[row, chrom]
  region = regions[row, region]
  region_start = regions[row, region_start]
  region_end = regions[row, region_end]
  region_pos = regions[row, pos]
  
  for(phenotype in phenotypes){
  
    data_tmp = phenotype_data[phenotype][[1]]
    
    #qc filtering
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
    id_conv_file = file=paste0(path_cojo_input, "cojo_id_converstion", region, "_", phenotype, "_pos_", region_pos, "__chr", region_chrom,".csv")
    
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
    
    fwrite(cojo_input, file=paste0(path_cojo_input, "cojo_input_region_", region, "_", phenotype, "_pos_", region_pos, "__chr", region_chrom,".ma"), row.names = F, quote = F, sep = " ")
    
    
    
  }
}

# Build and execute COJO command for each created input file --------------
message("\n--------------------------\n")
message("Performing COJO...\n")
input_files = list.files(paste0(projectpath, path_cojo_input))
i=0
for(input_file in input_files){
  chrom = str_split(input_file,"__")[[1]][2]
  chrom = str_remove(chrom, ".ma")
  chrom = str_remove(chrom, "chr")
  output_file = str_replace(input_file, "input", "output")
  output_file = str_replace(output_file, "\\.ma", "")
  
  gctaCommand = str_glue("{gctaCall} --bfile {path_plink}chr_{chrom} --chr {chrom} --cojo-slct --maf 0.01 --cojo-p 5e-8 --cojo-file {path_cojo_input}{input_file} --out {path_cojo_output}{output_file}")
  system(gctaCommand)
  
  i=i+1
  message("\n--------------------------\n")
  message(str_glue("Done {i} of {length(input_files)} files...\n"))
}
message("\n--------------------------\n")
message("Finished COJO select.\n")
