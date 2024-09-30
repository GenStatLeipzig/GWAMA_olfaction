# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2023-06-26
#
# Script Description: Make RA plots of the regions in the locus definition
#
#
# Notes:
# pipeline_name: 06c_RA_plots_LIFE_LD_franz_publication.R
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
n.cores = 40

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

#####
# 3. load data for SNPs
#####

pheno_files_list = list.files(path_data, ".gz")

preprocess = function(file) {
  dat = fread(paste0(path_data, file), nThread = 1)
  
#quality filter
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
  paste0(path_ra_plots, "RA_plots_defined_loci_LIFE_LD_publication.pdf"), #TODO check name 
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
  rs = regions[l, rsID ]
  
  #get locus data
  locus = phenotype_data[myPheno][[1]]
  locus = locus[(chrom == myChrom) &
                  (pos >= myStart) & (pos <= myEnd),
                c("markerID", "pos", "pFEM", "chrom")]
  
  colnames(locus) = c("snp", "position", "pvalue", "chr")
  
  locus[, snp := paste0("chr", snp)]
  matched = match(locus[, snp], changedIDs[, ID.formatted])
  locus[!is.na(matched), snp := changedIDs[matched[!is.na(matched)], ID.original]]
  
  #adapt id of leadsnp if necessary
  if(mySnp %in% changedIDs$ID.formatted){
    mySnp.checked = changedIDs[ID.formatted==mySnp]$ID.original
  } else {
    mySnp.checked = mySnp
  }
  
  
  input = createInputfilesRegAssocPlot1kg38Ensembl(
    locus,
    r_on_server = T,
    gene_validation_level = c("KNOWN"),
    path_ldreference = ldreference_fn,
    leadsnp = mySnp.checked
  )
  
  #plot
  mySubSize = 0.65
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
message("Finished RA Plots.\n")

dev.off()
