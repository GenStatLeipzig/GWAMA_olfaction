#' ---
#' title: "Analysis for credible sets of independent loci"
#' subtitle: "Credible Sets"
#' author: "Katrin Horn"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
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
#' 
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
#' 
#' 

# Script does the credible set analyisis for all phenotypes that are associated with a region
# i.e. the phenotypes where the top variant was observed
# for snps with multiple variants per region, the statistics are overwritten by the conditioned
# data.
# 
# only works on Aman due to unresolved library issues
# pipeline_name: 10a_credible_set_analyses_conditional.R

#TODO be aware that IDs in cojo might have changed, see the id_conversion files in cojo_input
rm(list=ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
library(gtx)
setwd(projectpath)

#qc filters
maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

#TODO enter regions where COJO conditional was applied (as string 'regionID_phenotype_subgroup')
# cojo_regions = c("5_cinnamon_all") 
cojo_regions = c()

#' # Preparations ####
#' ***
#' ##  ####

#' ## Load loci information
loci = fread(paste0(path_locus_definition, "locus_definition.csv"), dec=",") # TODO check locus definition file
path_cojo_cond = "07a_cojo/output_cond/" #TODO link to the COJO conditional results

# get a list of all the phenotypes to analyse
phenotypes = as.character(str_split(loci$phenotypes_in_region, " \\| ", simplify = T))
phenotypes = unique(phenotypes)
phenotypes = phenotypes[phenotypes!=""]
# phenotypes = phenotypes[str_starts(phenotypes, "cinnamon|pineapple_all")]

#' ## Load data of phenotypes for SNPs which can be done with unconditional statistics
phenotype_data = foreach(phenotype = phenotypes) %do% {
  filename = list.files(path_data, pattern = paste0("_", phenotype, ".*\\.gz"))
  if(length(filename)<1 || length(filename)>1){
    message(str_glue("File for {phenotype} not found or not unique"))
  } else {
    dat = fread(paste0(path_data, filename), nThread = 10)
    
    #qc filtering
    dat = dat[nWeightedMAF > maf_filter &
                nWeightedInfoScore > info_filter &
                I2 < i2_filter  &
                numberStudies >= n_studies_filter &
                numberLargeStudies >= n_large_studies_filter,]
    
    dat
  }
}
names(phenotype_data) = phenotypes

#' # Definition of credible set function
#' ***
CredibleSetFunction = function(data,filenm) { 
  
  dat = copy(data)
  
  #' ### Get Prior
  quant_0.975 = quantile(x = dat$beta, probs = 0.975, na.rm = T)
  quant_0.975   
  quant_0.025 = quantile(x = dat$beta, probs = 0.025, na.rm = T)
  quant_0.025  
  range_quant = quant_0.975-quant_0.025
  prior = range_quant/(2*1.96)
  print(prior)
  print(qnorm(p = 0.975, sd = prior) - quant_0.975)
  
  #' Should be close to 0 ...
  #' 
  #' ### Calculate Bayes factors ####
  #' Bayes factors, sum of Bayes factors and posterior probability (PP)
  dat[,abf.Wakefield:=gtx::abf.Wakefield(beta = beta, se = beta_se, priorsd = prior)]
  
  sum.abf.r1 = sum(dat[,abf.Wakefield], na.rm=T)
  sum.abf.r1
  dat[, PostProb:=abf.Wakefield/sum.abf.r1]
  summary(dat[,PostProb])
  ordering = order(dat[,PostProb], decreasing = TRUE)
  dat = dat[ordering,]
  dat[, SumProb:=cumsum(PostProb)]
  dat[, prior := prior]
  dat
  
  #' ### Summary ####
  myfile = paste0("08_credible_set_analysis/",filenm,".txt")
  write.table(dat,file=myfile,col.names = T, row.names = F, quote = F)
  
  return(dat)
}

#' # Prepare input data and calculate credible set (unconditional statistics) ####
#' ***
neededCol = c("Chr", "SNP", "bp", "refA", "freq", "beta", "beta_se", "p", "n", "info")  
 
uncond = foreach(line = 1:nrow(loci)) %do% {
  
  #get variables needed for data preparation
  myRegion = loci[line, region]
  region.start = loci[line, region_start]
  region.stop = loci[line, region_end]
  region.chrom = loci[line, chrom]
  
  #TODO switch to phenotypes in region column if CS for all phenotypes with 
  # gwsig hits is wanted instead of only the top phenotype
  
  phenotypes_in_region = loci[line, phenotype]
  
  # phenotypes_in_region = loci[line,phenotypes_in_region]
  # phenotypes_in_region = as.character(str_split(phenotypes_in_region, " \\| ", simplify = T))
  
  report_list = data.table()
  for(pheno in phenotypes_in_region){
    s = pheno
    cojo_id = paste0(myRegion, "_", pheno)
    
    locusData = phenotype_data[[pheno]]
    locusData = locusData[chrom == region.chrom & pos >= region.start & pos <= region.stop, ]
    
    
    if(cojo_id %in% cojo_regions){
      # procedure for regions with conditional statistics
      
      # collect files
      cojo_cond_files = list.files(path_cojo_cond, pattern = cojo_id, full.names = T)
      cojo_cond_files = subset(cojo_cond_files, str_ends(cojo_cond_files, ".cma.cojo"))
      
      for (cond_file in cojo_cond_files){
        
        #overwrite unconditional statistics with conditional statistics
        cond_data = fread(cond_file)
        
        setnames(cond_data, c("bC", "bC_se", "pC"), c("betaFEM", "seFEM", "pFEM"),
                 skip_absent = T)
        cond_data[,SNP := str_remove(SNP, "chr")]
        
        cond_locus_data = copy(locusData)
        cond_locus_data[, betaFEM := NULL]
        cond_locus_data[, seFEM := NULL]
        cond_locus_data[, pFEM := NULL]
        
        # exclude snps that have no conditional statistics
        cond_locus_data = merge(cond_locus_data, cond_data[,c("SNP", "betaFEM", "seFEM", "pFEM")],
                                by.x = "markerID", by.y = "SNP", all = F)
        
        # define credible set
        setnames(
          cond_locus_data,
          c(
            "chrom",
            "markerID",
            "pos",
            "ea",
            "nWeightedEAF",
            "betaFEM",
            "seFEM",
            "pFEM",
            "totalN",
            "nWeightedInfoScore"
          ),
          neededCol
        )
        colsOut = setdiff(colnames(cond_locus_data), neededCol)
        cond_locus_data[, get("colsOut") := NULL]
        
        # snp = loci[line, markerID]
        snp = str_match(cond_file, "chr(.*)\\.cma\\.cojo$")[2]
        result = CredibleSetFunction(cond_locus_data, paste0("region_",myRegion, "_",s, "_", gsub(":", "_", snp)))
        
        report = list(
          phenotype = pheno,
          region = myRegion,
          SNP = snp,
          CredSet.95 = nrow(result[SumProb <= 0.95,]) + 1,
          CredSet.99 = nrow(result[SumProb <= 0.99,]) + 1,
          prior = unique(result[, prior])
        )
        report = as.data.table(report)
        report_list = rbind(report_list, report)
        
      }
      
    } else {
      # procedure for unconditional statistics
      
      setnames(locusData, c("chrom", "markerID", "pos", "ea", "nWeightedEAF", "betaFEM", "seFEM", "pFEM", "totalN", "nWeightedInfoScore"), neededCol)
      colsOut = setdiff(colnames(locusData), neededCol)
      locusData[, get("colsOut") := NULL]
      
      snp = loci[line, markerID]
      result = CredibleSetFunction(locusData, paste0("region_",myRegion, "_",s, "_", gsub(":", "_", snp)))
      
      report = list(phenotype = pheno, region = myRegion, SNP = snp, CredSet.95 = nrow(result[SumProb <= 0.95, ])+1, CredSet.99 = nrow(result[SumProb <= 0.99, ])+1, prior = unique(result[, prior]))
      report = as.data.table(report)
      report_list = rbind(report_list, report)
    }
  }
  report_list
}

uncond = rbindlist(uncond)

#' Save results in file
overview = uncond
myOrder = order(overview[, phenotype], overview[, region])
overview = overview[myOrder, ]

#get summary of prior used
summary(overview[, prior])

write.table(overview, file = "08_credible_set_analysis/SNP_numbers_Credible_Sets.txt", col.names = T, row.names = F, quote = F)

#' # Generate files for annotation with GWAS Pipeline ####
#' ***
colsNeeded = c("SNP", "CredSet", "PostProb", "SumProb")
phenotypes = unique(overview[, phenotype])

lists = foreach(p = phenotypes) %do% {
  myFiles = list.files(path = "08_credible_set_analysis/", pattern = paste0("_", p))

  dat = foreach(f = myFiles) %do% {
    myFile = fread(paste0("08_credible_set_analysis/",f))
    pos = min(which(myFile[,SumProb]>0.99))
    myFile = myFile[c(1:pos), ]
    dummy = unlist(strsplit(f, split = "_"))[2]
    myFile[, CredSet := paste0("Region", dummy)]
    
    colsOut = setdiff(colnames(myFile), colsNeeded)
    myFile[, get("colsOut"):=NULL]
    setnames(myFile, "SNP", "snp")
    setcolorder(myFile, c("snp", "CredSet", "PostProb", "SumProb"))
    myFile
  }
  dat = rbindlist(dat)
  outFile = paste0("08_credible_set_analysis/", p, "_SNPs_to_annotate.txt")
  write.table(dat, file = outFile, col.names = T, row.names = F, quote = F, sep = "\t")
}

lists = rbindlist(lists)


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

