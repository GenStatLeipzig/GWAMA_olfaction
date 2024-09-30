#' ---
#' title: "Mendelian Randomization 4 - FDR and plots"
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
#' This is the forth MR script for MR4. Here, I will load the results, and perform FDR by exposure and setting (to discuss with Markus and Franz). 
#' 
#' Once I have the significant exposure - outcome pairs, I will create the following plots: 
#' 
#' - scatter plot per combination
#' - forest plots of all combinations at once
#' 
#' In addition, I repeat the MR analyses using other methods to check how robust the causal association is. 
#'  
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR4_04_getFDR_plots.R

rm(list = ls())
time0<-Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(paste0(projectpath, "00_scripts/"))
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
load(paste0(path_MR_results,MR_prefix,"/03_MRresults_IVW.RData"))

#' # Loop for FDR ####
#' ***
#' The disease data was not sex-stratified. I use "all" as main analysis, and the stratified data as sensitivity checks. 
#' 
MRTab = copy(MRTab_MR4)
MRTab[,setting := gsub(".*_","",outcome)]
MRTab[,outcome := gsub("_.*","",outcome)]
mySettings = unique(MRTab$setting)

dumTab1 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  MRTab2 = copy(MRTab)
  MRTab2 = MRTab2[setting==mySettings[i]]
  
  # create a heatmap of the unadjusted pvalues
  dumTab <- dcast(MRTab2, exposure ~ outcome, value.var="pval_IVW")
  dumMat = as.matrix(dumTab[,2:14])
  rownames(dumMat) = dumTab$exposure
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  
  fn = paste0(path_MR_results,MR_prefix,"/04_PvalueHeatmap_unadjusted_",mySettings[i],".png")
  png(filename = fn,width = 2000, height = 2000, res=250)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # perform multiple testing correction
  myFDR = addHierarchFDR(pvalues = MRTab2[,pval_IVW], 
                         categs = MRTab2[,exposure],
                         quiet = F)
  MRTab2[,pval_adjusted := myFDR$fdr_level1]
  MRTab2[,hierFDR := myFDR$hierarch_fdr5proz]
  
  # create a heatmap of the adjusted pvalues
  dumTab <- dcast(MRTab2, exposure ~ outcome, value.var="pval_adjusted")
  dumMat = as.matrix(dumTab[,2:14])
  rownames(dumMat) = dumTab$exposure
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  
  fn = paste0(path_MR_results,MR_prefix,"/04_PvalueHeatmap_adjusted_",mySettings[i],".png")
  png(filename = fn,width = 2000, height = 2000, res=250)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  MRTab2
}

MRTab_FDR = rbindlist(dumTab1)

MRTab_sig = copy(MRTab_FDR)
MRTab_sig = MRTab_sig[pval_adjusted<0.05,]

MRTab_nominal = copy(MRTab_FDR)
MRTab_nominal = MRTab_nominal[pval_IVW<0.05,]

#' **Summary**
#' 
#' 16 nominal significant results, 4 still significant after FDR correction: 
#' 
#' - ALZHEIMER on SCORE in all
#' - AD_WIDE on lemon in women
#' - ALZHEIMER on lemon in women
#' - AD_WIDE on liquorice in men
#'  
#' # Forest Plots ####
#' ***
#' I want forest plots for ALZHEIMER and AD_WIDE on SCORE, lemon and liquorice in all sex-settings
#' 
plotData = copy(MRTab_FDR)
plotData = plotData[outcome %in% c("SCORE","lemon","liquorice"),]
plotData = plotData[exposure %in% c("ALZHEIMER","AD_WIDE"),]

myOutcomes = unique(plotData$outcome)

dumTab2 = foreach(i=1:length(myOutcomes))%do%{
  #i=3
  plotData1 = copy(plotData)
  plotData1 = plotData1[outcome == myOutcomes[i],]
  setorder(plotData1,beta_IVW)
  
  plotData1[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  plotData1[,upperCI95 := beta_IVW+1.96*SE_IVW]
  xmin = min(c(0,plotData1$lowerCI95, plotData1$upperCI95), na.rm = T)
  xmax = max(c(0,plotData1$lowerCI95, plotData1$upperCI95), na.rm = T)
  xmin2 = -0.5
  xmax2 = 0.5
  if(i==3) xmin2 = -0.15
  if(i==3) xmax2 = 0.05
  xrange = xmax2 - xmin2
  xmin_margin = xmin2 - 0.7*xrange
  xmax_margin = xmax2 + 0.7*xrange
  
  myHeader = paste0("Risk of Alzheimers disease on ",myOutcomes[i]," recognition") 
  if(i==3)  myHeader = paste0("Risk of Alzheimers disease on smelling ",myOutcomes[i]) 
  myRows = paste0("in ",plotData1$setting," - ",plotData1$exposure," (",plotData1[,NR_SNPs_total]," SNPs)")
  
  fn = paste0(path_MR_results,MR_prefix,"/04_ForestPlot_Alzheimers_",myOutcomes[i],".png")
  png(filename = fn,width = 2400, height = 1200, res=250)
  par(mar=c(5,6,0,4))
  par(font=1)
  dets = forest(x=plotData1[,beta_IVW], 
                sei = plotData1[,SE_IVW],
                header = myHeader,
                slab = myRows,
                ylim = c(-0, nrow(plotData1)+3) , 
                cex =1, 
                xlim = c(xmin_margin, xmax_margin), 
                alim = c(xmin2, xmax2))
  dev.off()
  
}

#' # Scatter plots and different methods ####
#' ***
#' For Alzheimers and the three outcomes, I want the scatter plots and other MR methods. 
myExposures = unique(plotData$exposure)

ToDoList = data.table(outcome = rep(myOutcomes,each=2),
                      exposure = rep(myExposures,3),
                      setting = rep(c("female","male","all"),each=2))

load(paste0(path_MR_results,MR_prefix,"/01_SNP_on_neuroDiseases.RData"))
load(paste0(path_MR_results,MR_prefix,"/02_SNP_on_olfactoryPerception.RData"))

#' check outcome data: I want all SNPs to be available in all outcomes! There is one variant that was filtered in the sex-stratified data due to MAC problems. 
#' 
test = myTab_OP[,.N, by=rsID]
test[N!=39,]
myTab_filtered[rsID == test[N!=39,rsID]]
myTab_OP = myTab_OP[!is.element(rsID,test[N!=39,rsID]),]
myTab_filtered = myTab_filtered[!is.element(rsID,test[N!=39,rsID]),]

dumTab3 = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=6
  myRow = ToDoList[i,]
  
  data_Y = copy(myTab_OP)
  data_Y = data_Y[phenotype == myRow$outcome & setting == myRow$setting,]
  
  data_X = copy(myTab_filtered)
  data_X = data_X[disease == myRow$exposure,]
  
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
                    exposure = myRow$exposure,
                    outcome = paste0(myRow$outcome," in ",myRow$setting),
                    snps = paste(data_X$rsID,data_X$cyto,sep="\n"))
  
  myPlot = mr_plot(mr_obj, interactive = F,labels = T)
  fn = paste0(path_MR_results,MR_prefix,"/04_ScatterPlot_",myRow$exposure,"_",myRow$outcome,"_",myRow$setting,".png")
  png(filename = fn,width = 2000, height = 2000, res=250)
  print(myPlot)
  dev.off()
  
  fn = paste0(path_MR_results,MR_prefix,"/04_LeaveOneOut_",myRow$exposure,"_",myRow$outcome,"_",myRow$setting,".png")
  png(filename = fn,width = 2000, height = 2000, res=250)
  mr_loo(mr_obj)
  dev.off()
  
  res2 = mr_allmethods(mr_obj)
  res3 = res2@Values
  setDT(res3)
  names(res3) = c("Method","Estimate","SE","lowerCI95","upperCI95","pvalue")
  res3[,exposure := myRow$exposure]
  res3[,outcome := myRow$outcome]
  res3[,setting_sex := myRow$setting]
  res3[,NR_SNPs_total := dim(data_X)[1]]
  res3
}

MR_otherMethods = rbindlist(dumTab3)

#' # Save results ####
#' ***

WriteXLS(x = c("MRTab_FDR","MR_otherMethods"), 
         ExcelFileName=paste0(path_MR_results,MR_prefix,"/04_MR_NeuroDiseases_Smelling.xlsx"), 
         SheetNames=c("IVW_allCombis","otherMethods_sigOnly"), 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(MRTab_FDR,file = paste0(path_MR_results,MR_prefix,"/04_MRresults_IVW_FDRadj.RData"))
save(MR_otherMethods,file = paste0(path_MR_results,MR_prefix,"/04_MRresults_otherMethods.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
