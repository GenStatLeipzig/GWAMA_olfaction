#' ---
#' title: "Mendelian Randomization 2 - run MR analyses with multiple instruments"
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
#' This is the third script for MR2. Here, I will run the MR analyses per exposure (different sniffin sticks) and each outcome (3 neuro diseases). 
#' 
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR2_03_MR_IVW.R

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

#' # Load data ####
#' ***
load(paste0(path_MR_results,MR_prefix,"/02_MRbySNP_Smelling_NeuroDiseases.RData"))
myMRresults2[,dumID := paste0(exposure,"_", exposure_setting,"__", outcome)]
ToDoList = myMRresults2[,.N,dumID]
ToDoList = ToDoList[N>1,] 

myTab1 = fread(paste0("../", path_locus_definition, "/locus_definition_rsID.csv"), dec = ",")
myNames = c("rsID","markerID","chrom","pos","aa","ea","nWeightedMAF","nWeightedInfoScore",
            "phenotype","betaFEM","seFEM","pFEM","I2","totalN")
stopifnot(myNames %in% names(myTab1))
colsOut<-setdiff(colnames(myTab1),myNames)
myTab1[,get("colsOut"):=NULL]
setcolorder(myTab1,myNames)
myTab1[,gene := c("GSX2/FIP1L1","ADCY2","FBXL17","OR2","TAAR5","OR51/OR52","","OR5","OR5","","OR11")]

matched = match(myMRresults2$rsID,myTab1$rsID)
myMRresults2[,gene:= myTab1[matched,gene]]
myMRresults2[,dumID2 := paste0(rsID,"\n",gene)]

#' # Loop #### 
#' ***
#' 
dumTab1 = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  myData = copy(myMRresults2)
  myData = myData[dumID == myRow$dumID,]
  
  # create MR object 
  mr_obj = mr_input(bx = myData$exposure_beta,
                    bxse = myData$exposure_SE,
                    by = myData$outcome_beta,
                    byse = myData$outcome_SE,
                    exposure = gsub("__.*","",myRow$dumID),
                    outcome = gsub(".*__","",myRow$dumID),
                    snps = myData$dumID2)
  res1 = mr_ivw(mr_obj)
  res2 = data.table(exposure = res1@Exposure,
                    outcome = res1@Outcome,
                    NR_SNPs_total = res1@SNPs,
                    beta_IVW = res1@Estimate,
                    SE_IVW = res1@StdError,
                    pval_IVW = res1@Pvalue,
                    Fstat = res1@Fstat,
                    HeteroStat = res1@Heter.Stat[1],
                    HeteroStat_pval = res1@Heter.Stat[2],
                    comment = "all SNPs")
  if(res2$pval_IVW<0.05){
    myPlot = mr_plot(mr_obj, interactive = F,labels = T)
    fn = paste0(path_MR_results,MR_prefix,"/03_ScatterPlot_all_",res1@Exposure,"__",res1@Outcome,".png")
    png(filename = fn,width = 2000, height = 2000, res=250)
    print(myPlot)
    dev.off()
    
    fn = paste0(path_MR_results,MR_prefix,"/03_LeaveOneOut_all_",res1@Exposure,"__",res1@Outcome,".png")
    png(filename = fn,width = 2000, height = 2000, res=250)
    mr_loo(mr_obj)
    dev.off()
    
  }
  
  res2
}
myMRresults3 = rbindlist(dumTab1)
myMRresults3[,table(pval_IVW<0.05)]
myMRresults3[pval_IVW<0.1]
myMRresults3[,pval_adj := p.adjust(pval_IVW,method="fdr"),by = c("exposure","comment")]
myMRresults3[pval_adj<0.1]
myMRresults3[pval_adj == min(pval_adj) & pval_IVW == min(pval_IVW),]

#' I will create scatter plots for the best IVW result: cinnamon (all) on ALZHEIMERS
#' 
myData = copy(myMRresults2)
myData = myData[dumID == "cinnamon_female__ALZHEIMER",]

# create MR object 
mr_obj = mr_input(bx = myData$exposure_beta,
                  bxse = myData$exposure_SE,
                  by = myData$outcome_beta,
                  byse = myData$outcome_SE,
                  exposure = "cinnamon_female",
                  outcome = "ALZHEIMER",
                  snps = myData$dumID2)
res1 = mr_ivw(mr_obj)
myPlot = mr_plot(mr_obj, interactive = F,labels = T)
fn = paste0(path_MR_results,MR_prefix,"/03_ScatterPlot_female_",res1@Exposure,"__",res1@Outcome,".png")
png(filename = fn,width = 2000, height = 2000, res=250)
print(myPlot)
dev.off()
  
fn = paste0(path_MR_results,MR_prefix,"/03_LeaveOneOut_female_",res1@Exposure,"__",res1@Outcome,".png")
png(filename = fn,width = 2000, height = 2000, res=250)
mr_loo(mr_obj)
dev.off()

#' # Save results ####
#' ***

WriteXLS(x = c("myMRresults3"), 
         ExcelFileName=paste0(path_MR_results,MR_prefix,"/03_MRIVW_Smelling_NeuroDiseases.xlsx"), 
         SheetNames=c("IVW"), 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(myMRresults3,file = paste0(path_MR_results,MR_prefix,"/03_MRIVW_Smelling_NeuroDiseases.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
