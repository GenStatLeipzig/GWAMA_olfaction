#' ---
#' title: 'Colocalisation Analysis: Male versus Female'
#' author: 'Andreas Kühnapfel, Franz Förster'
#' date: 'Last compiled on `r format(Sys.time(), '%d %B %Y')`'
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' Test if SNPs colocalise between the sexes. 
#'
#'
#' # Initialize ####
#' ***
#' Update server & set working directory
# pipeline_name: 11b_coloc_sex_ia.R

#TODO change type when analysing the score

rm(list = setdiff(ls(), "allinfo"))
time0 = Sys.time()


source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
setwd(projectpath)

library(coloc)

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

#TODO select regionIDs without significant sex effect 
# regions_to_analyze = c()  

#####################################################
########## PERFORM COLOCALISATION ANALYSIS ##########
#####################################################

##### Setup
#

##### Common Base Data
locus = fread(
  paste0(
    path_locus_definition,
    "locus_definition.csv" #TODO check locus definition file
  ),
  dec = ","
)
# locus = locus[which(locus$region %in% regions_to_analyze),] # TODO uncomment if only specific regions are wanted

data = data.frame(
  snp = locus$markerID,
  chr = locus$chrom,
  pos = locus$pos,
  locus_ll = locus$region_start,
  locus_ul = locus$region_end,
  region = locus$region,
  trait = locus$phenotype
)
data$trait = str_replace(data$trait, "(_all|_female|_male)", "")

data_1_lifted = data
data_2_lifted = data

coloc = list()
for (i in 1:nrow(locus)) {
  phenotype = locus[i, phenotype]
  phenotype = str_replace(phenotype, "(_all|_female|_male)", "")
  
  # load summary data corresponding to phenotype
  meta_file_all = list.files(paste0(path_data), pattern = paste0("GWASMA_", phenotype, "_all.*\\.gz"))
  meta_file_male = list.files(paste0(path_data), paste0("GWASMA_", phenotype, "_male.*\\.gz"))
  meta_file_female = list.files(paste0(path_data), pattern = paste0("GWASMA_", phenotype, "_female.*\\.gz"))
  
  # sumStat_1     = fread(paste0(path_data, meta_file_all), nThread = 30 ) # all results not needed
  sumStat_2 = fread(paste0(path_data, meta_file_male), nThread = 30)
  sumStat_3 = fread(paste0(path_data, meta_file_female), nThread = 30)
  
  ### Create Locus-Wise Lists as Input for Colocalisation Analyses
  
  ### Step/Option 1: Use step10 Data
  
  step10_1_lifted = sumStat_2  # select male data
  
  filt_1_locus_i = (step10_1_lifted$chrom == data_1_lifted$chr[i]) &
    (step10_1_lifted$pos >= data_1_lifted$locus_ll[i]) &
    (step10_1_lifted$pos <= data_1_lifted$locus_ul[i]) &
    (step10_1_lifted$numberStudies >= n_studies_filter) &
    (step10_1_lifted$numberLargeStudies >= n_large_studies_filter) & 
    (step10_1_lifted$I2 < i2_filter) &
    (step10_1_lifted$nWeightedInfoScore > info_filter) &
    (step10_1_lifted$nWeightedMAF > maf_filter)


  
  step10_1_lifted_filt_1_locus_i = step10_1_lifted[filt_1_locus_i,]
  
  #check for homogeneous sample size
  # sample size is now inhomogenous due to variing sample sizes in ARIC
  # stopifnot(length(unique(step10_1_lifted_filt_1_locus_i$totalN)) == 1)  #Sample size should be the same for all SNPs
  
  gwas_1_i = list(
    pvalues = step10_1_lifted_filt_1_locus_i$pFEM,
    totalN = as.numeric(step10_1_lifted_filt_1_locus_i$totalN),
    MAF = step10_1_lifted_filt_1_locus_i$nWeightedMAF,
    beta = step10_1_lifted_filt_1_locus_i$betaFEM,
    varbeta = (step10_1_lifted_filt_1_locus_i$seFEM) ^ 2,
    type = rep("cc", (sum(filt_1_locus_i, na.rm = TRUE) + sum(is.na(filt_1_locus_i)))),
    snp = step10_1_lifted_filt_1_locus_i$markerID
  )
  
  ### Step/Option 1: Use step10 Data ->
  step10_2_lifted = sumStat_3  # select female data
  
  filt_2_locus_i = (step10_2_lifted$chrom == data_2_lifted$chr[i]) &
    (step10_2_lifted$pos >= data_2_lifted$locus_ll[i]) &
    (step10_2_lifted$pos <= data_2_lifted$locus_ul[i]) &
    (step10_2_lifted$numberStudies >= n_studies_filter) &
    (step10_2_lifted$numberLargeStudies >= n_large_studies_filter) &
    (step10_2_lifted$I2 < i2_filter) & 
    (step10_2_lifted$nWeightedMAF > maf_filter)
    (step10_2_lifted$nWeightedInfoScore > info_filter)
  
  step10_2_lifted_filt_2_locus_i = step10_2_lifted[filt_2_locus_i,]
  
  gwas_2_i = list(
    pvalues = step10_2_lifted_filt_2_locus_i$pFEM,
    totalN = as.numeric(step10_2_lifted_filt_2_locus_i$totalN),
    MAF = step10_2_lifted_filt_2_locus_i$nWeightedMAF,
    beta = step10_2_lifted_filt_2_locus_i$betaFEM,
    varbeta = (step10_2_lifted_filt_2_locus_i$seFEM) ^ 2,
    type = rep("cc", (sum(filt_2_locus_i, na.rm = TRUE) + sum(is.na(filt_2_locus_i)))),
    snp = step10_2_lifted_filt_2_locus_i$markerID
  )
  
  ### Perform Colocalisation Analyses
  
  coloc_ij = list()
  
  dum_k_1 = is.na(gwas_1_i$pvalues) |
    is.na(gwas_1_i$totalN) | 
    is.na(gwas_1_i$MAF) |
    is.na(gwas_1_i$beta) |
    is.na(gwas_1_i$varbeta) |
    is.na(gwas_1_i$type) | is.na(gwas_1_i$snp)
  dum_k_2 = is.na(gwas_2_i$pvalues) |
    is.na(gwas_2_i$totalN) |
    is.na(gwas_2_i$MAF) |
    is.na(gwas_2_i$beta) |
    is.na(gwas_2_i$varbeta) |
    is.na(gwas_2_i$type) | is.na(gwas_2_i$snp)
  
  gwas_1_filtered = gwas_1_i
  gwas_2_filtered = gwas_2_i
  for (l in 1:length(gwas_1_filtered)) {
    
    gwas_1_filtered[[l]] = gwas_1_filtered[[l]][!dum_k_1]
    gwas_2_filtered[[l]] = gwas_2_filtered[[l]][!dum_k_2]
  }
  
  dum_k = match(gwas_1_filtered[["snp"]], gwas_2_filtered[["snp"]])
  for (l in 1:length(gwas_1_filtered)) {
    
    gwas_1_filtered[[l]] = gwas_1_filtered[[l]][!is.na(dum_k)]
    gwas_2_filtered[[l]] = gwas_2_filtered[[l]][dum_k[!is.na(dum_k)]]
  }
  
  #set type and N to a single value as the coloc function expects
  gwas_1_filtered[["type"]] = gwas_1_filtered[["type"]][1]
  gwas_1_filtered[["totalN"]] = gwas_1_filtered[["totalN"]][1]
  gwas_2_filtered[["type"]] = gwas_2_filtered[["type"]][1]
  gwas_2_filtered[["totalN"]] = gwas_2_filtered[["totalN"]][1]
  
  
  # check_dataset(gwas_1_filtered)
  # check_dataset(gwas_2_filtered)
  coloc_ij = coloc.abf(gwas_1_filtered, gwas_2_filtered)
  
  
  
  ### Save Colocalisation Analysis Results
  coloc[[i]] = coloc_ij
}

##### Write Results
res_summary = data.frame()

for (i in 1:length(coloc)) {
  res_summary = rbind(res_summary, coloc[[i]]$summary)
}
names(res_summary) = c("nsnps",
                       "PP.H0.abf",
                       "PP.H1.abf",
                       "PP.H2.abf",
                       "PP.H3.abf",
                       "PP.H4.abf")
res_summary = data.frame(
  locus = paste0("Region ", data$region),
  trait1 = paste0(data$trait, "_male"),
  trait2 = paste0(data$trait, "_female"),
  res_summary
)

write.table(
  res_summary,
  file = "11_sex_interaction/coloc_sex_ia.csv",
  sep = ";",
  dec = ",",
  row.names = FALSE
)


########## END OF FILE ##########

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "mins"), 3),
        " minutes")

