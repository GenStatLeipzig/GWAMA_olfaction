# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster and Katrin Horn
# 
# Date: 2024-01-11
#
# Script Description: Make RA plots of region 7 by manually specifying all necessary pairs of SNPs for LD visualisation
# 
#
# Notes:
# pipeline_name: 06d_RA_plot_LIFE_LD_Region7_cleaned_publication.R
#

#####
# 0.
#####
rm(list = ls())

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
require(ggplot2) #plotting
require(plyr) #data manipulation
require(reshape) #data manipulation
require(stringr) #string manipulation
require(Hmisc) #frank harrell's miscellaneous functions
ht <- function ( d, myrows=10 ) rbind ( head ( d , myrows ), tail ( d , myrows ))
setwd(projectpath)
n.cores = 10

#####
# 1. define plot function and important variables (function from Holger (GWAS Pipeline))
#####
source("helper_scripts/RegAssocPlot_hg38_EnsemblGenes_LIFE_LD.R")

locus_definition_file = "locus_definition_rsID.csv" #TODO check file


maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

#####
# 2. load region specifications from locus definition
#####

regions = fread(paste0(path_locus_definition, locus_definition_file), dec = ",")
regions = regions[7,]
regions

#####
# 3. load data for SNPs
#####

pheno_files_list = list.files(path_data, ".gz")
pheno_files_list = pheno_files_list[31] 
pheno_files_list

preprocess = function(file) {
  dat = fread(paste0(path_data, file), nThread = 1)
  
#quality filter
#TODO comment in when filtering is wanted
  dat = dat[nWeightedMAF > maf_filter &
              nWeightedInfoScore > info_filter &
              I2 < i2_filter  &
              numberStudies >= n_studies_filter &
              numberLargeStudies >= n_large_studies_filter,]
  return(dat)
}
my.cluster <- parallel::makeCluster(n.cores,
                                    type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

#parallel processing
message("\n--------------------------\n")
message("loading SNP data...\n")
phenotype_data = foreach (file = pheno_files_list,
                          .packages = c("data.table", "stringr")) %dopar% {
                            preprocess(file)
                          }
parallel::stopCluster(cl = my.cluster)

getPheno = function(f) {
  phenotype = str_match(f, "GWASMA_(\\D+)_\\d")[2]
  return(phenotype)
}

names(phenotype_data) = lapply(pheno_files_list, getPheno)

#####
# 4. plot function with data generation routine
#####
.env = new.env() #to prevent an error in the plotting script that expects an environment
myRows = c(1:nrow(regions))

pdf(
  paste0(path_ra_plots, "RA_plot_region_7_aa_excluded_special.pdf"),
  width = 14,
  height = 7
)

changedIDs = fread("additional_information/changed_life_ids.txt")



done = foreach(l = myRows) %do% {
  
  myPheno = regions[l, phenotype]
  myPmin = regions[l, pFEM]
  mySnp = paste0("chr", regions[l, markerID])
  myRegion = regions[l, region]
  myChrom = regions[l, chrom]
  myStart = regions[l, region_start]
  myEnd = regions[l, region_end]
  myMAF = round(regions[l, nWeightedMAF], 2)
  myPos = regions[l, pos]
  myWidth = ((myEnd - myStart)/2) %/% 1000
  myInfo = round(regions[l, nWeightedInfoScore],2)
  rs = regions[l, rsID]
  
  #get locus data
  locus = phenotype_data[myPheno][[1]]
  locus = locus[(chrom == myChrom) &
                  (pos >= myStart) & (pos <= myEnd),
                c("markerID", "pos", "pFEM", "chrom")]
  
  colnames(locus) = c("snp", "position", "pvalue", "chr")
  
  locus[, snp := paste0("chr", snp)]
  matched = match(locus[, snp], changedIDs[, ID.formatted])
  locus[!is.na(matched), snp := changedIDs[matched[!is.na(matched)], ID.original]]
  
  input = createInputfilesRegAssocPlot1kg38Ensembl(
    locus,
    r_on_server = T,
    gene_validation_level = c("KNOWN"),
    path_ldreference = ldreference_fn,
    leadsnp = mySnp
  )
  
  myLocus = input$myLocus
  
  #calculate LD for every SNP with leadSNP
  mySNPs = myLocus$NAME
  length(mySNPs)
  ld = foreach(s = mySNPs) %do% {
    sname = gsub(pattern = ":", replacement = "_", s)
    plinkCall = paste0(plink2, " --pfile ", pgenFile, " --ld ", input$leadsnp, " ", s, " --out 00_scripts/tmp/06d_LD_result_", sname)
    # system(plinkCall) #TODO uncomment when calculation is wanted
    
    resFile = paste0("00_scripts/tmp/06d_LD_result_", sname, ".log")
    sysCall = paste0("grep r^2 ", resFile)
    res = system(sysCall, intern = T)  
    if (length(grepl(pattern = "r^2", x = res)) > 0) {
      dummy = unlist(strsplit(res, split = "  "))[2]
      rsquare = as.numeric(unlist(strsplit(dummy, split = "="))[2])
    } else {
      rsquare = NA
    }
    
    result = list(SNP = s, R2 = rsquare)
    as.list(result)
  }
  ld = rbindlist(ld)
  
  #change the LD (Rsquare data in input)
  matched = match(input$myLocus$NAME, ld[, SNP])
  sum(is.na(matched))   #should be 0
  input$myLocus$RSQR = ld[matched, R2]
  
  #plot
  mySubSize = 0.65
  # myRegion = 7 #TODO after final decision of validity of loci
  customASplot_woBROAD(
    locus = input$myLocus,
    lead_snp = input$leadsnp,
    map = input$myMap,
    genes = input$myGenes,
    shownregion_kb = myWidth,
    maintitle = paste0("Region ", myRegion , ": ", rs),
    subtitle = paste0("trait: ", myPheno, ";  rsq: ", myInfo, ";  MAF: ", myMAF),
    weakR2 = 0.2,
    center_lead_snp = F,
    cex_genname = 0.5,
    title_size = 1,
    subtitle_size = mySubSize,
    gene_lines = 5,
    sf = c(2, 3),
    logpmax = -log10(myPmin)
  )
  
  as.list(mySnp)
}
done = rbindlist(done)

message("\n--------------------------\n")
message("Finished RA Plot.\n")

dev.off()
