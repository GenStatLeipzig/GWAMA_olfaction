#' ---
#' title: "Mendelian Randomization 1 - extract SNP effects on OP"
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
#' This is the second script for MR1. Here, I will get the SNP - outcome association (SH --> OP). The SNP list has been generated in the previous script.
#'  
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR1_02_SNP_Outcome.R

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

#' # Load data ####
#' ***
#' ## Load SNP list 
load(paste0(path_MR_results,MR_prefix,"/01_SNPList_instruments.RData"))
myTab[,chrPosEAOA := paste(CHR,POS,EA,OA,sep=":")]
myTab[,chrPosOAEA := paste(CHR,POS,OA,EA,sep=":")]

#' ## Load outcome data 
statistics = list.files(path = path_data, pattern = ".gz")
statistics2 = gsub("GWASMA_","",statistics)
statistics2 = gsub("_2024.*","",statistics2)
statistics3 = unlist(strsplit(statistics2,"_"))

ToDoList = data.table(phenotype = statistics3[seq(1,length(statistics3),2)],
                      setting = statistics3[seq(2,length(statistics3),2)],
                      filename = statistics)

myNames_old = c("phenotype","setting",
                "markerID","chrom","pos","ea","aa","nWeightedEAF","nWeightedMAF","nWeightedInfoScore","I2",
                "betaFEM","seFEM","pFEM","zscore","totalN")
myNames_new = c("phenotype","setting",
                "markerID","CHR","POS","EA","OA","EAF","MAF","info","I2",
                "beta","SE","pvalue","zscore","sampleSize")

dumTab1 = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myTime = Sys.time()
  myRow = ToDoList[i,]
  message("\nWorking on phenotype ",myRow$phenotype," in ",myRow$setting," (",i," of ",dim(ToDoList)[1],")")
  
  # load data
  myfn1 = paste0(path_data,myRow$filename)
  erg1 = fread(myfn1)
  erg1[,phenotype := myRow$phenotype]
  erg1[,setting := myRow$setting]
  erg1[,zscore := betaFEM/seFEM]
  
  # filter columns
  stopifnot(myNames_old %in% names(erg1))
  colsOut<-setdiff(colnames(erg1),myNames_old)
  erg1[,get("colsOut"):=NULL]
  setcolorder(erg1,myNames_old)
  names(erg1) = myNames_new

  # filter rows
  erg1 = erg1[markerID %in% c(myTab$chrPosOAEA,myTab$chrPosEAOA)]
  erg1[,chrPos := paste(CHR,POS,sep=":")]
  setorder(erg1,CHR,POS)
  
  myTab2 = copy(myTab)
  myTab2[,chrPos := paste(CHR,POS,sep=":")]
  matched = match(erg1$chrPos,myTab2$chrPos)
  myTab2 = myTab2[matched,]
  stopifnot(sum(duplicated(myTab2$chrPos))==0)
  
  # harmonize alleles
  stopifnot(erg1$CHR == myTab2$CHR)
  stopifnot(erg1$POS == myTab2$POS)
  erg1[,rsID := myTab2$rsID,]
  filt = myTab2$EA != erg1$EA
  erg1[filt,EAF := 1-EAF]
  erg1[filt,beta := (-1)*beta]
  erg1[,diff := abs(EAF - myTab2$EAF)]
  erg1 = erg1[diff<0.2]
  
  # get rid of some matching and filtering columns
  erg1[,diff:= NULL]
  
  # save in temp file
  save(erg1,file = paste0(path_MR_temp,"MR1_",myRow$phenotype,"_",myRow$setting,"_filtered.RData"))
  erg1
}
myTab_OP = rbindlist(dumTab1)
save(myTab_OP,file = paste0(path_MR_results,MR_prefix,"/02_SNP_on_olfactoryPerception.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
