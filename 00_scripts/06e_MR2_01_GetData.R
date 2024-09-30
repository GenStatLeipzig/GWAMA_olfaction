#' ---
#' title: "Mendelian Randomization 2 - load and match data"
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
#' This is the first script for MR2. Here, I will load all genome-wide significant SNPs (Table 1) and extract their statistics from the score trait. 
#'   
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR2_01_GetData.R

rm(list = ls())
time0<-Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(paste0(projectpath, "00_scripts/"))
.libPaths()

MR_prefix = "MR2_Smelling_on_neuroDiseases"
if(dir.exists(paste0(path_MR_results,MR_prefix))==F){
  dir.create(paste0(path_MR_results,MR_prefix))
  message("Created results folder ",paste0(path_MR_results,MR_prefix))
}else{
  message("Using pre-existing results folder ",paste0(path_MR_results,MR_prefix))
}

#' # Get Instruments ####
#' ***
myTab1 = fread(paste0("../",path_locus_definition,"/locus_definition_rsID.csv"), dec = ",")
myNames = c("rsID","markerID","chrom","pos","aa","ea","nWeightedMAF","nWeightedInfoScore",
            "phenotype","betaFEM","seFEM","pFEM","I2","totalN")
stopifnot(myNames %in% names(myTab1))
colsOut<-setdiff(colnames(myTab1),myNames)
myTab1[,get("colsOut"):=NULL]
setcolorder(myTab1,myNames)
myTab1

#' Feedback from Franz: rs11245786 should be excluded as it is not completely independent of rs669453 (LD R2 = 0.1941 in EUR)
myTab1 = myTab1[!grepl("rs11245786",rsID),]

allPhenos = unique(myTab1$phenotype)
allPhenos = gsub("_.*","",allPhenos)
allPhenos = unique(allPhenos)

statistics = list.files(path = path_data, pattern = ".gz")
statistics2 = gsub("GWASMA_","",statistics)
statistics2 = gsub("_2024.*","",statistics2)
statistics3 = unlist(strsplit(statistics2,"_"))

ToDoList = data.table(phenotype = statistics3[seq(1,length(statistics3),2)],
                      setting = statistics3[seq(2,length(statistics3),2)],
                      filename = statistics)
ToDoList = ToDoList[phenotype %in% allPhenos,]

myNames_old = c("phenotype","setting",
                "markerID","chrom","pos","ea","aa","nWeightedEAF","nWeightedMAF","nWeightedInfoScore","I2",
                "betaFEM","seFEM","pFEM","zscore","totalN")
myNames_new = c("phenotype","setting",
                "markerID","CHR","POS","EA","OA","EAF","MAF","info","I2",
                "beta","SE","pvalue","zscore","sampleSize")

dumTab1 = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  message("\nWorking on phenotype ",myRow$phenotype," in ",myRow$setting," (",i," of ",dim(ToDoList)[1],")")
  
  # load data
  myfn1 = paste0(path_data,myRow$filename)
  erg1 = fread(myfn1)
  erg1 = erg1[markerID %in% myTab1$markerID,]
  matched = match(myTab1$markerID,erg1$markerID)
  erg1 = erg1[matched,]
  erg1[,phenotype := myRow$phenotype]
  erg1[,setting := myRow$setting]
  erg1[,zscore := betaFEM/seFEM]
  
  # filter columns
  stopifnot(myNames_old %in% names(erg1))
  colsOut<-setdiff(colnames(erg1),myNames_old)
  erg1[,get("colsOut"):=NULL]
  setcolorder(erg1,myNames_old)
  names(erg1) = myNames_new
  
  stopifnot(erg1$EA == myTab1$ea)
  stopifnot(erg1$OA == myTab1$aa)
  
  erg1
}
myTab2 = rbindlist(dumTab1)
save(myTab2,file = paste0(path_MR_results,MR_prefix,"/01_SNP_on_olfactoryPerception.RData"))

myTab2[,table(pvalue<1e-6,setting)]
myTab2[,table(pvalue<0.05,setting)]

#' Add EAF from SCORE all to myTab1 (which is used for harmonization)
matched = match(myTab1$markerID,myTab2$markerID)
myTab1[,EAF := myTab2[matched,EAF]]

#' # Load neuro outcome #####
#' ***
statistics = list.files(path = path_NeuroData, pattern = ".gz")
statistics
neuroDiseases = gsub("ukbb_summary_stats_G6_","",statistics)
neuroDiseases = gsub("_meta_out.tsv.gz","",neuroDiseases)

myTab1[,chrPosEAOA := paste(chrom,pos,ea,aa,sep=":")]
myTab1[,chrPosOAEA := paste(chrom,pos,aa,ea,sep=":")]

dumTab2 = foreach(i = 1:3)%do%{
  #i=1
  
  # load and filter data
  myfn2 = paste0(path_NeuroData,statistics[i])
  erg2 = fread(myfn2)
  erg2 = erg2[rsid %in% myTab1$rsID,]
  
  stopifnot(myTab1$chrom == erg2$`#CHR`)
  stopifnot(myTab1$pos == erg2$POS)
  
  # check effect allele frequency
  plot(myTab1$EAF,erg2$FINNGEN_af_alt)
  abline(0,1)
  plot(myTab1$EAF,erg2$UKBB_af_alt)
  abline(0,1)
  
  # check minor allele frequency
  erg2[,FINNGEN_maf_alt := FINNGEN_af_alt]
  erg2[,UKBB_maf_alt := UKBB_af_alt]
  erg2[FINNGEN_maf_alt>0.5,FINNGEN_maf_alt := 1-FINNGEN_maf_alt]
  erg2[UKBB_maf_alt>0.5,UKBB_maf_alt := 1-UKBB_maf_alt]
  plot(myTab1$nWeightedMAF,erg2$FINNGEN_maf_alt)
  abline(0,1)
  plot(myTab1$nWeightedMAF,erg2$UKBB_maf_alt)
  abline(0,1)
  
  # check alleles
  table(myTab1$ea == erg2$ALT)
  table(myTab1$aa == erg2$ALT)
  filt = myTab1$aa==erg2$ALT
  
  # harmonize alleles 
  erg2[, beta := all_inv_var_meta_beta]
  erg2[filt, beta := beta * (-1)]
  erg2[, SE := all_inv_var_meta_sebeta]
  erg2[, pval := all_inv_var_meta_p]
  erg2[, logP := all_inv_var_meta_mlogp]
  erg2[, heteroStat := all_inv_var_het_p]
  erg2[, EAF_FG := FINNGEN_af_alt]
  erg2[filt, EAF_FG := 1-EAF_FG]
  erg2[, EAF_UKBB := UKBB_af_alt]
  erg2[filt, EAF_UKBB := 1-EAF_UKBB]
  erg2 = erg2[,c(22,5,1:4,25:31)]
  
  # check effect allele frequency again
  plot(myTab1$EAF,erg2$EAF_FG)
  abline(0,1)
  plot(myTab1$EAF,erg2$EAF_UKBB)
  abline(0,1)
  
  # return data
  erg2[,outcome := neuroDiseases[i]]
  erg2[,ALT := myTab1[,ea]]
  erg2[,REF := myTab1[,aa]]
  erg2[,matchingID := myTab1$markerID]
  erg2
}
myTab3 = rbindlist(dumTab2)

setnames(myTab3,"rsid","rsID")
setnames(myTab3,"#CHR","CHR")
setnames(myTab3,"REF","OA")
setnames(myTab3,"ALT","EA")

save(myTab3,file = paste0(path_MR_results,MR_prefix,"/01_SNP_on_NeuroOutcomes.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
