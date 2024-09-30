# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2023-06-22
#
# Script Description: Source File that is loaded in all scripts to provide common paths
#
# Notes:
# pipeline_name: 00_SourceFile_smelling_meta.R
# 
#

Sys.umask(mode = "0002") #make files readable for all

#############################
# common paths
#############################
projectpath = "GWAMA_olfaction/" #TODO set project directory (should be ...GWAMA_olfaction/)
path_raw_data = "01_metal_output/" # Summary statistics in GWAS cataloge SSF format
path_data = paste0(projectpath, "02_MetaGWAS/") # summary statistics formatted for use in pipeline 
path_mh_plots = "04_mh_qq_plots_overview_stats/" 
path_locus_definition = "05_locus_definition/" 
path_ra_plots = "06_ra_plots/"
path_snp_annotation = "" #TODO ensembl annotation path, data taken from https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/ (downloaded at 20.07.2023)
gctaCall = "gcta_1.92.0beta3/gcta64" #TODO link gcta package location for COJO 
path_plink = "" #TODO plink files for usage in COJO (LIFE-Adult was used as reference)
ldreference_fn = "" #TODO link LD reference for RA plots

path_chainfile = "hg38ToHg19.over.chain" #TODO link chain file hg38tohg19

p_gtx = "2007_GTEx_v8/" # TODO path to GTEx data

plink2 = "plink2.0/20240105/plink2" #TODO link to plink2 executable
pgenFile = "" #TODO link to pgen file (needed for recalculation of LD based on LIFE-Adult with plink for redoing RA plot of region 7)

# Path for MR analyses
path_MR_results = "../12_MR/"
path_MR_temp = "../12_MR_temp/"
path_NeuroData = "../additional_information/2023_FinnGen_UKBB_NeuroPhenotypes/"
path_SHdata = "../additional_information/2020_Ruth_SexHormones_UKBB/"
path_KEGGgenes = "../additional_information/KEGG_SteroidHormoneBiosynthesis/"


#############################
# R library and R packages that are commonly used across all scripts
#############################
.libPaths("") #TODO set libpath depending on your computer
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(formatR))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(MendelianRandomization))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(toolboxH))
suppressPackageStartupMessages(library(corrplot))

#############################
# helper functions
#############################

calculate_maf_from_freq <- function(x) {
  if (x < 0.5) {
    return(x)
  } else {
    return(1 - x)
  }
}

getSmallestDist = function(x) {
  if(length(x)>1){
    y = c(x[2:length(x)], max(x)+1000000)
    z = min(y-x)
  }else{
    z=1000000
  }
  return(z)
}

posBasedPruning = function(data, distance=1e6,colname_chr="CHR",colname_pos="POS",colname_logP="logP"){
  
  stopifnot(colname_chr %in% names(data))
  stopifnot(colname_pos %in% names(data))
  stopifnot(colname_logP %in% names(data))
  
  if(colname_chr != "CHR") setnames(data,colname_chr,"CHR")
  if(colname_pos != "POS") setnames(data,colname_pos,"POS")
  if(colname_logP != "logP") setnames(data,colname_logP,"logP")
  
  myCHRs = unique(data$CHR)
  
  result.22 = foreach(s2 = myCHRs) %do% {
    # s2 = myCHRs[1]
    subdata2 = copy(data)
    subdata2 = subdata2[CHR == s2, ]
    setkey(subdata2, POS)
    
    if(dim(subdata2)[1]<=1){
      subdata2[, keep := T]
      subdata2[, NR_SNPs := 0]
    }else{
      subdata2[, keep := NA]
      subdata2[, NR_SNPs := as.numeric(NA)]
      
      smallestDist = getSmallestDist(subdata2[, POS])
      while(smallestDist < distance) {
        #minP = min(subdata2[is.na(keep), p])
        maxLogP = max(subdata2[is.na(keep), logP])
        myPOS = subdata2[maxLogP == logP & is.na(keep), POS]
        if(length(myPOS)>1){
          myPOS = myPOS[1]
        }
        subdata2[POS == myPOS, keep := T]
        
        #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
        myFilt = (subdata2[, POS] < (myPOS - distance)) | 
          (subdata2[, POS] > (myPOS + distance)) | 
          subdata2[, keep] 
        myFilt[is.na(myFilt)] = FALSE
        subdata2 = subdata2[myFilt == TRUE, ]
        
        subdata2[POS == myPOS, NR_SNPs := sum(myFilt==F)]
        smallestDist = getSmallestDist(subdata2[, POS])
      }
      
      #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
      subdata2[is.na(keep), NR_SNPs := 0]
      subdata2[is.na(keep), keep := TRUE]
    }
    
    subdata2
  }
  data_filtered = rbindlist(result.22)
  setorder(data_filtered,CHR,POS)
  
  return(data_filtered)
}


MRfunction_jp<-function(betaX,seX,betaY,seY){
  betaIV<-betaY/betaX
  se.ratio.st<- seY/sqrt(betaX^2)
  se.ratio.nd<- sqrt(seY^2/betaX^2 + betaY^2*seX^2/betaX^4)
  p1<-2*pnorm(-abs(betaIV/se.ratio.st))
  p2<-2*pnorm(-abs(betaIV/se.ratio.nd))
  res<-data.table(
    beta_IV=betaIV,
    se_IV1=se.ratio.st,
    se_IV2=se.ratio.nd,
    p_IV1=p1,
    p_IV2=p2)
  return(res)
}
