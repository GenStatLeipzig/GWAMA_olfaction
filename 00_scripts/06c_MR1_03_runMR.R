#' ---
#' title: "Mendelian Randomization 1 - run MR analyses"
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
#' This is the third script for MR1. Here, I will run the MR analyses per exposure (different sex hormones) and each outcome (39 OPs). 
#' 
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR1_03_runMR.R

rm(list = ls())
time0<-Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(paste0(projectpath, "00_scripts/"))
.libPaths()

MR_prefix = "MR1_steroidHormone_on_Smelling"
if(dir.exists(paste0(path_MR_results,MR_prefix))==F){
  dir.create(paste0(path_MR_results,MR_prefix))
  message("Created results folder ",paste0(path_MR_results,MR_prefix))
}else{
  message("Using pre-existing results folder ",paste0(path_MR_results,MR_prefix))
}

#' # Loop: all pathway genes ####
#' ***
load(paste0(path_MR_results,MR_prefix,"/01_SNP_on_steroidHormones.RData"))
load(paste0(path_MR_results,MR_prefix,"/02_SNP_on_olfactoryPerception.RData"))

#' Filter for genome-wide significant in at least one setting (in the first script, I only filtered for suggestive significance)
test = myTab_filtered[,max(abs(zscore)),by=rsID]
test = test[V1<5.45,]
myTab_filtered = myTab_filtered[!is.element(rsID,test$rsID)]

myTab_filtered[, exposure := paste(hormone,setting,sep="_")]
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
    res4[,SNPset := "sameSNPs"]
    
    # rerun with significant SNPs (Z-score > 5) - different SNP sets for men - women - combined!
    if(myExposures[i]!="E2_Men"){
      filt = abs(data_X$zscore)>5.45
      mr_obj = mr_input(bx = data_X$beta[filt],
                        bxse = data_X$SE[filt],
                        by = data_Y$beta[filt],
                        byse = data_Y$SE[filt],
                        exposure = myExposures[i],
                        outcome = myOutcomes[j])
      res2 = mr_ivw(mr_obj)
      res5 = data.table(exposure = myExposures[i],
                        outcome = myOutcomes[j],
                        NR_SNPs_total = res2@SNPs,
                        beta_IVW = res2@Estimate,
                        SE_IVW = res2@StdError,
                        pval_IVW = res2@Pvalue,
                        FStat = res2@Fstat,
                        HeteroStat = res2@Heter.Stat[1],
                        HeteroStat_pval = res2@Heter.Stat[2])
      res5[,SNPset := "strongSNPs"]
      
      # rerun with shared significant SNPs
      filt = grepl("shared",data_X$flag)
      mr_obj = mr_input(bx = data_X$beta[filt],
                        bxse = data_X$SE[filt],
                        by = data_Y$beta[filt],
                        byse = data_Y$SE[filt],
                        exposure = myExposures[i],
                        outcome = myOutcomes[j])
      res2 = mr_ivw(mr_obj)
      res7 = data.table(exposure = myExposures[i],
                        outcome = myOutcomes[j],
                        NR_SNPs_total = res2@SNPs,
                        beta_IVW = res2@Estimate,
                        SE_IVW = res2@StdError,
                        pval_IVW = res2@Pvalue,
                        FStat = res2@Fstat,
                        HeteroStat = res2@Heter.Stat[1],
                        HeteroStat_pval = res2@Heter.Stat[2])
      res7[,SNPset := "sameStrongSNPs"]
      
      res6 = rbind(res4,res5,res7)
      
    }else{
      res6 = res4
    }
    res6
  }
  tab2 = rbindlist(dumTab2)
  tab2
}
MRTab_Pathway = rbindlist(dumTab1)
MRTab_Pathway[pval_IVW<0.05,]

save(MRTab_Pathway,file = paste0(path_MR_results,MR_prefix,"/03_MRresults_IVW.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
