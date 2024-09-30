# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2023-06-26
#
# Script Description: Make RA plots of the regions in the locus definition. Plots not only the RA Plot of the top phenotype but also the two other subgroups (all, male, female). Marks the SNPs within the region with a coloured circle.
# Uses LD calculated with PLINK for the plots instead of TOPLD
#
#
# Notes:
# pipeline_name: 06e_RA_plots_triplets_with_snp_comparison_LIFE_LD.R
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
regions = regions[c(6),] #TODO select region

#TODO restrict widow

# restrictions for region 6
regions[1, region_start:= 5400000]
regions[1, region_end := 5550000]


#####
# 3. load data for SNPs
#####

pheno_files_list = list.files(path_data, ".gz")
# pheno_files_list = pheno_files_list[1] 

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

changedIDs = fread("additional_information/changed_life_ids.txt") # list of snps where IDs changed during study harmonization to match them to the LIFE reference

pdf(
  paste0(path_ra_plots, "RA_plots_regions_with_idependant_variants.pdf"), #TODO check name
  width = 8,
  height = 8
)
layout(matrix(c(1,1,2,3), nrow=2, byrow = T))


done = foreach(l = myRows) %do% {
  
  myPheno = regions[l, phenotype]
  myPmin = regions[l, pFEM]
  # mySnp = regions[l, markerID]
  myRegion = regions[l, region]
  myChrom = regions[l, chrom]
  myStart = regions[l, region_start]
  myEnd = regions[l, region_end]
  # myMAF = round(regions[l, nWeightedMAF], 2)
  # myPos = regions[l, pos]
  myWidth = ((myEnd - myStart)/2) %/% 1000
  # myInfo = round(regions[l, nWeightedInfoScore],2)
  
  trait = str_split(myPheno, "_", simplify = T)[1]
  myPhenos = paste0(trait, c("_all", "_male", "_female"))
  
  #collect the top snps for each subgroup
  
  top_snps = foreach (myPheno = myPhenos) %do% {
    locus = phenotype_data[myPheno][[1]]
    locus = locus[(chrom == myChrom) &
                    (pos >= myStart) & (pos <= myEnd),]
    
    #collect snp data
    snp_data = locus[pFEM == min(locus$pFEM),]
    snp_data = snp_data[1,] #prevents errors when multiple snps have same p-value
    mySnp = snp_data$markerID
  }
  top_snp_a = paste0("chr",top_snps[[1]])
  top_snp_m = paste0("chr",top_snps[[2]])
  top_snp_f = paste0("chr", top_snps[[3]])
  
  #adapt IDs if necessary
  if(top_snp_a %in% changedIDs$ID.formatted){
    top_snp_a = changedIDs[ID.formatted==top_snp_a]$ID.original
  }
  if(top_snp_m %in% changedIDs$ID.formatted){
    top_snp_m = changedIDs[ID.formatted==top_snp_m]$ID.original
  }
  if(top_snp_f %in% changedIDs$ID.formatted){
    top_snp_f = changedIDs[ID.formatted==top_snp_f]$ID.original
  }
  
  
  for(myPheno in myPhenos){
    #get locus data
    locus = phenotype_data[myPheno][[1]]
    locus = locus[(chrom == myChrom) &
                    (pos >= myStart) & (pos <= myEnd),]
    
    #collect snp data
    snp_data = locus[pFEM == min(locus$pFEM),]
    snp_data = snp_data[1,] #prevents errors when multiple snps have same p-value
    mySnp = paste0("chr", snp_data$markerID)
    myMAF = snp_data$nWeightedMAF
    myPos = snp_data$pos
    myInfo = round(snp_data$nWeightedInfoScore, 2)
    
    locus = locus[,c("markerID", "pos", "pFEM", "chrom")]
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
    
    # do not plot reference genes
    myGenes = input$myGenes
    myGenes = myGenes[0,]
    
    #some specific formatting for the publication
    if (str_detect(myPheno, "_all")){
      maintitle = paste0("Regional Association Plot of locus ", myRegion, " - ", str_replace(myPheno, "_", " "))
    } else {
      maintitle = str_replace(myPheno, "_", " ")
    }
    
    # TODO
    # rsIDs for locus 6
    snp_to_rs = c("chr11:5504313:T:C" = "rs317787", "chr11:5471572:G:A" = "rs430332")
    # rsIDs for locus 11
    # snp_to_rs = c("chr14:20201081:C:T" = "rs2318888", "chr14:20217151:T:C" = "rs8007085", "chr14:20198925:A:G" = "rs11159353")
    
    
    
    #plot
    mySubSize = 0.65
    customASplot_woBROAD(
      locus = input$myLocus,
      lead_snp = input$leadsnp,
      map = input$myMap,
      genes = myGenes,
      shownregion_kb = myWidth,
      maintitle = maintitle,
      subtitle = paste0("SNP: ", snp_to_rs[mySnp.checked], ";  rsq: ", myInfo, ";  MAF: ", myMAF),
      weakR2 = 0.2,
      center_lead_snp = F,
      cex_genname = 0.5,
      title_size = 1,
      subtitle_size = mySubSize,
      gene_lines = 1,
      sf = c(10, 9),
      logpmax = -log10(myPmin)
    )
    
    #add points for the different variants
    points(x=locus[snp == top_snp_a, position], y= -log10(locus[snp == top_snp_a, pvalue]),
           pch = 21, col = "black", cex=3, lw = 3) # not used for region 6 because it is identical to the top-SNP
    points(x=locus[snp == top_snp_m, position], y= -log10(locus[snp == top_snp_m, pvalue]),
           pch = 21, col = "green", cex =3, lw = 3)
    points(x=locus[snp == top_snp_f, position], y= -log10(locus[snp == top_snp_f, pvalue]),
           pch = 21, col = "magenta", cex = 3, lw = 3)
  }
  as.list(mySnp)
}
done = rbindlist(done)

message("\n--------------------------\n")
message("Finished RA Plots.\n")

dev.off()
