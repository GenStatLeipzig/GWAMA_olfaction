#' ---
#' title: "Mendelian Randomization 4 - extract instruments for AD and PD"
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
#' This is the first script for MR4. Here, I will get the SNP - exposure association (AD/PD --> OP). 
#'    
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR4_01_GetInstruments.R

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

#' # Load neuro outcome #####
#' ***
statistics = list.files(path = path_NeuroData, pattern = ".gz")
statistics
neuroDiseases = gsub("ukbb_summary_stats_G6_","",statistics)
neuroDiseases = gsub("_meta_out.tsv.gz","",neuroDiseases)
ToDoList = data.table(pheno = neuroDiseases,
                      stats = statistics)

#' Sample sizes obtained from https://public-metaresults-fg-ukbb.finngen.fi/, filtered for Alzheimers and Parkinsons
ToDoList[, N_cases := c(15617 + 839, 10520 + 831, 4681 + 1660)]
ToDoList[, N_control := c(396564 + 410833, 401661 + 419700, 407500 + 403249)]
ToDoList[, N_total := c(15617 + 839 + 396564 + 410833, 10520 + 831 + 401661 + 419700, 4681 + 1660 + 407500 + 403249)]

erg1 = fread(paste0(path_data,"GWASMA_pineapple_all_2024-03-01.gz"))

myNames_old = c("disease","setting",
                "rsid","chrPosOAEA","chrPosEAOA","CHR","POS","ALT","REF","FINNGEN_af_alt","UKBB_af_alt",
                "all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p","zscore","N_total","N_cases","N_controls")
myNames_new = c("disease","setting",
                "rsID","chrPosOAEA","chrPosEAOA","CHR","POS","EA","OA","EAF_FinnGen","EAF_UKBB",
                "beta","SE","pvalue","zscore","N_total","N_cases","N_controls")

dumTab2 = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  message("Working on disease ",ToDoList[i,pheno])
  
  # load and filter data
  myfn2 = paste0(path_NeuroData,ToDoList[i,stats])
  erg2 = fread(myfn2)
  
  # I only want SNPs which were available in both UKB and FINNGEN
  erg2 = erg2[all_meta_N==2,]
  
  # I only want genome-wide significant SNPs with low heterogeneity
  erg2 = erg2[all_inv_var_meta_p<=5e-8,]
  erg2 = erg2[all_inv_var_het_p > 0.05,]
  
  # check if SNPs are in OP GWAS data
  setnames(erg2,"#CHR","CHR")
  erg2[,chrPosOAEA := paste(CHR, POS, REF, ALT, sep=":")]
  erg2[,chrPosEAOA := paste(CHR, POS, ALT, REF, sep=":")]
  filt = is.element(erg2$chrPosOAEA,erg1$markerID) | is.element(erg2$chrPosEAOA,erg1$markerID)
  erg2 = erg2[filt,]
  
  # priority pruning by position filtering
  erg3 = posBasedPruning(data = erg2,distance = 1000000,colname_chr = "CHR",
                         colname_pos = "POS",colname_logP = "all_inv_var_meta_mlogp")
  
  # select relevant columns
  erg3[,disease := ToDoList[i,pheno]]
  erg3[,setting := "sex-combined"]
  erg3[,N_cases := ToDoList[i,N_cases]]
  erg3[,N_controls := ToDoList[i,N_control]]
  erg3[,N_total := ToDoList[i,N_total]]
  erg3[, zscore := all_inv_var_meta_beta / all_inv_var_meta_sebeta]
  
  stopifnot(myNames_old %in% names(erg3))
  colsOut<-setdiff(colnames(erg3),myNames_old)
  erg3[,get("colsOut"):=NULL]
  setcolorder(erg3,myNames_old)
  names(erg3) = myNames_new
  erg3
}
myTab_filtered = rbindlist(dumTab2)

save(myTab_filtered,file = paste0(path_MR_results,MR_prefix,"/01_SNP_on_neuroDiseases.RData"))

#' # SNP List - updated ####
#' ***
#' In the next script, I have to load and filter all 12 x 3 OP data sets. I need a SNPList to do so. 
#' 
myTab = copy(myTab_filtered)
myTab = myTab[!duplicated(rsID),]
myTab = myTab[,c(3:11)]
myTab[,CHR2 := CHR]
myTab[,CHR := as.numeric(CHR)]
myTab[CHR2=="X",CHR := 23]
setorder(myTab,CHR,POS)

save(myTab,file = paste0(path_MR_results,MR_prefix,"/01_SNPList_instruments.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
