#' ---
#' title: "Mendelian Randomization 4 - run MR analyses"
#' subtitle: "GWAS olfactory perception"
#' author: "Janne Pott"
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
#' This is the third script for MR4. Here, I will run the MR analyses per exposure (different neurodegenerative diseases) and each outcome (39 OPs). 
#' 
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR4_03_runIVW.R

rm(list = ls())
time0<-Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(paste0(projectpath, "00_scripts"))
.libPaths()

MR_prefix = "MR3_neuroDiseases_on_Smelling"
if(dir.exists(paste0(path_MR_results,MR_prefix))==F){
  dir.create(paste0(path_MR_results,MR_prefix))
  message("Created results folder ",paste0(path_MR_results,MR_prefix))
}else{
  message("Using pre-existing results folder ",paste0(path_MR_results,MR_prefix))
}

#' # Load data ####
#' ***
load(paste0(path_MR_results,MR_prefix,"/01_SNP_on_neuroDiseases.RData"))
load(paste0(path_MR_results,MR_prefix,"/02_SNP_on_olfactoryPerception.RData"))

#' check outcome data: I want all SNPs to be available in all outcomes! There is one variant that was filtered in the sex-stratified data due to MAC problems. 
#' 
test = myTab_OP[,.N, by=rsID]
test[N!=39,]
myTab_filtered[rsID == test[N!=39,rsID]]
myTab_OP = myTab_OP[!is.element(rsID,test[N!=39,rsID]),]
myTab_filtered = myTab_filtered[!is.element(rsID,test[N!=39,rsID]),]

#' add some filtering here
#' 
myTab_filtered[,absZ := abs(zscore)]
myTab_filtered[,Fstat := zscore^2]
myTab_filtered[, exposure := disease]
myExposures = unique(myTab_filtered$exposure)
myTab_OP[, outcome := paste(phenotype,setting,sep="_")]
myOutcomes = unique(myTab_OP$outcome)

dumTab1 = foreach(i = 1:length(myExposures))%do%{
  #i=1
  data_X = copy(myTab_filtered)
  data_X = data_X[exposure == myExposures[i],]
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    #j=1
    data_Y = copy(myTab_OP)
    data_Y = data_Y[outcome == myOutcomes[j]]
    
    # match data
    matched = match(data_X$rsID,data_Y$rsID)
    data_Y = data_Y[matched,]
    filt = is.na(data_Y$beta) | is.na(data_X$beta)
    data_Y = data_Y[!filt,]
    data_X = data_X[!filt,]
    stopifnot(data_X$rsID == data_Y$rsID)
    
    # create MR object 
    mr_obj = mr_input(bx = data_X$beta,
                      bxse = data_X$SE,
                      by = data_Y$beta,
                      byse = data_Y$SE,
                      exposure = myExposures[i],
                      outcome = myOutcomes[j])
    res2 = mr_ivw(mr_obj)
    res4 = data.table(exposure = myExposures[i],
                      outcome = myOutcomes[j],
                      NR_SNPs_total = res2@SNPs,
                      beta_IVW = res2@Estimate,
                      SE_IVW = res2@StdError,
                      pval_IVW = res2@Pvalue,
                      FStat = res2@Fstat,
                      HeteroStat = res2@Heter.Stat[1],
                      HeteroStat_pval = res2@Heter.Stat[2])
    res4[,SNPset := "allSNPs"]
    res4
  }
  tab2 = rbindlist(dumTab2)
  tab2
}
MRTab_MR4 = rbindlist(dumTab1)
MRTab_MR4[pval_IVW<0.05,]

save(MRTab_MR4,file = paste0(path_MR_results,MR_prefix,"/03_MRresults_IVW.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
