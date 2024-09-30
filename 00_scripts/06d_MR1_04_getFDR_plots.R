#' ---
#' title: "Mendelian Randomization 1 - FDR and plots"
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
#' This is the forth MR script for MR1. Here, I will load the results, and perform FDR by exposure and setting. 
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
#' pipeline_name: 15_MR1_04_getFDR_plots.R

rm(list = ls())
time0<-Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(paste0(projectpath, "00_scripts"))
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
load(paste0(path_MR_results,MR_prefix,"/03_MRresults_IVW.RData"))

MRTab = copy(MRTab_Pathway)

#' Filter for the right combinations: all on all, women on women, men on men only
MRTab = MRTab[grepl("Comb",exposure) & grepl("all",outcome) | 
                grepl("_Women",exposure) & grepl("female",outcome) | 
                grepl("_Men",exposure) & grepl("_male",outcome)]


#' # Loop for FDR ####
#' ***
mySettings = unique(MRTab$SNPset)

dumTab1 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  MRTab2 = copy(MRTab)
  MRTab2 = MRTab2[SNPset==mySettings[i]]
  
  # create a heatmap of the unadjusted pvalues
  MRTab2[,outcome := gsub("_.*","",outcome)]
  dumTab <- dcast(MRTab2, exposure ~ outcome, value.var="pval_IVW")
  dummy = dumTab$exposure
  dummy = gsub("SHBG_BMIadj","SHBGadj",dummy)
  dummy = unlist(strsplit(dummy,"_"))
  dummy1 = dummy[seq(2,length(dummy),2)]
  dummy2 = dummy[seq(1,length(dummy),2)]
  dumTab$exposure2 = paste(dummy1,dummy2,sep=" - ")
  ordering = order(dumTab$exposure2)
  dumTab = dumTab[ordering,]
  dumMat = as.matrix(dumTab[,2:13])
  rownames(dumMat) = dumTab$exposure2
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  
  # fn = paste0(path_MR_results,MR_prefix,"/04_PvalueHeatmap_unadjusted_",mySettings[i],".png")
  # png(filename = fn,width = 2000, height = 2000, res=250)
  # corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  # dev.off()
  
  # perform multiple testing correction
  myFDR = addHierarchFDR(pvalues = MRTab2[,pval_IVW], 
                         categs = MRTab2[,exposure],
                         quiet = F)
  MRTab2[,pval_adjusted := myFDR$fdr_level1]
  MRTab2[,hierFDR := myFDR$hierarch_fdr5proz]
  
  # create a heatmap of the adjusted pvalues
  dumTab <- dcast(MRTab2, exposure ~ outcome, value.var="pval_adjusted")
  dummy = dumTab$exposure
  dummy = gsub("SHBG_BMIadj","SHBGadj",dummy)
  dummy = unlist(strsplit(dummy,"_"))
  dummy1 = dummy[seq(2,length(dummy),2)]
  dummy2 = dummy[seq(1,length(dummy),2)]
  dumTab$exposure2 = paste(dummy1,dummy2,sep=" - ")
  ordering = order(dumTab$exposure2)
  dumTab = dumTab[ordering,]
  dumMat = as.matrix(dumTab[,2:13])
  rownames(dumMat) = dumTab$exposure2
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  
  # fn = paste0(path_MR_results,MR_prefix,"/04_PvalueHeatmap_adjusted_",mySettings[i],".png")
  # png(filename = fn,width = 2000, height = 2000, res=250)
  # corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  # dev.off()
  
  MRTab2
}

MRTab_FDR = rbindlist(dumTab1)

MRTab_sig = copy(MRTab_FDR)
MRTab_sig = MRTab_sig[pval_adjusted<0.05,]
MRTab_sig

MRTab_nominal = copy(MRTab_FDR)
MRTab_nominal = MRTab_nominal[pval_IVW<0.05,]
MRTab_nominal[SNPset == "strongSNPs",]

#' **Summary**
#' 
#' There are 9 nominal significant estimates in the main analysis (strong SNPs): 
#' 
#' - BAT on fish in men
#' - BAT on pineapple in men
#' - SHBG on licorice in men
#' - SHBG on SCORE in men
#' - SHBG (BMI adj) on coffee in sex-combined setting
#' - SHBG on coffee in sex-combined setting
#' - TT on pineapple in men
#' - TT on banana in sex-combined setting
#' - TT on coffee in sex-combined setting
#'
#' Significant after FDR correction: 
#' 
#' - **BAT on pineapple in men**
#'  
#' # Forest Plots ####
#' ***
#' I want a forest plot for BAT and TT on pineapple including all sexes.
#' add FStat in Forest Plot?
plotData = copy(MRTab_FDR)
plotData = plotData[outcome == "pineapple",]
plotData = plotData[grepl("BAT",exposure) | grepl("TT",exposure),]

mySettings = unique(plotData$SNPset)
mySettings

dumTab2 = foreach(i=1:length(mySettings))%do%{
  #i=2
  plotData1 = copy(plotData)
  plotData1 = plotData1[SNPset == mySettings[i],]
  setorder(plotData1,beta_IVW)
  
  plotData1[,lowerCI95 := beta_IVW-1.96*SE_IVW]
  plotData1[,upperCI95 := beta_IVW+1.96*SE_IVW]
  xmin = min(c(0,plotData1$lowerCI95, plotData1$upperCI95), na.rm = T)
  xmax = max(c(0,plotData1$lowerCI95, plotData1$upperCI95), na.rm = T)
  xmin2 = floor(xmin)
  xmax2 = ceiling(xmax)
  if(xmax2>5) xmax2 = 5
  if(xmin2< (-5)) xmin2 = -5
  xrange = xmax2 - xmin2
  xmin_margin = xmin2 - 0.7*xrange
  xmax_margin = xmax2 + 0.7*xrange
  
  if(i==1) mySet = "Sensitivity 1: union of all instruments"
  if(i==2) mySet = "Main analysis: strong instruments per setting"
  if(i==3) mySet = "Sensitivity 2: colocalizing loci with strong instruments"
  myHeader = paste0("(Bioavailable) Testosterone levels on pineapple recognition") 
  myLabels = paste0(gsub("_"," - ",plotData1$exposure)," \n",plotData1[,NR_SNPs_total]," SNPs; F-Stat: ",plotData1[,round(FStat,1)])
  
  fn = paste0(path_MR_results,MR_prefix,"/04_ForestPlot_Testo_Pineapple_",mySettings[i],".png")
  png(filename = fn,width = 2400, height = 1500, res=250)
  par(mar=c(5,6,0,4))
  par(font=1)
  dets = forest(x=plotData1[,beta_IVW], 
                sei = plotData1[,SE_IVW],
                header = myHeader,
                slab = myLabels,
                ylim = c(-0, nrow(plotData1)+3) , 
                cex =1, 
                xlim = c(xmin_margin, xmax_margin), 
                alim = c(xmin2, xmax2))
  par(font=4)
  text(min(dets$xlim), max(dets$ylim-1.5), mySet, pos=4)
  dev.off()
  
}

#' # Scatter plots and different methods ####
#' ***
#' For BAT in men, I want the scatter plots and other MR methods. 

load(paste0(path_MR_results,MR_prefix,"/01_SNP_on_steroidHormones.RData"))
load(paste0(path_MR_results,MR_prefix,"/02_SNP_on_olfactoryPerception.RData"))

mySettings

dumTab3 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  mySetting = mySettings[i]
  
  data_Y = copy(myTab_OP)
  data_Y = data_Y[phenotype == "pineapple" & setting =="male",]
  
  data_X = copy(myTab_filtered)
  data_X = data_X[hormone == "BAT" & setting == "Men",]
  
  matched = match(data_X$rsID,data_Y$rsID)
  data_Y = data_Y[matched,]
  filt = is.na(data_Y$beta) | is.na(data_X$beta)
  data_Y = data_Y[!filt,]
  data_X = data_X[!filt,]
  
  if(mySetting == "strongSNPs"){
    filt = abs(data_X$zscore)>5.45
    table(filt)
    data_X = data_X[filt,]
    data_Y = data_Y[filt,]
  }else if(mySetting == "sameStrongSNPs"){
    filt = grepl("shared",data_X$flag)
    table(filt)
    data_X = data_X[filt,]
    data_Y = data_Y[filt,]
  }
  stopifnot(data_X$rsID == data_Y$rsID)
  
  # create MR object 
  mr_obj = mr_input(bx = data_X$beta,
                    bxse = data_X$SE,
                    by = data_Y$beta,
                    byse = data_Y$SE,
                    exposure = "BAT in men",
                    outcome = "pineapple in men",
                    snps = paste(data_X$rsID,data_X$cyto,sep="\n"))
  
  myPlot = mr_plot(mr_obj, interactive = F,labels = T)
  fn = paste0(path_MR_results,MR_prefix,"/04_ScatterPlot_BAT_pineapple_men_",mySetting,".png")
  png(filename = fn,width = 2000, height = 2000, res=250)
  print(myPlot)
  dev.off()
  
  fn = paste0(path_MR_results,MR_prefix,"/04_LeaveOneOut_BAT_pineapple_men_",mySetting,".png")
  png(filename = fn,width = 2000, height = 2000, res=250)
  mr_loo(mr_obj)
  dev.off()
  
  res2 = mr_allmethods(mr_obj)
  res3 = res2@Values
  setDT(res3)
  names(res3) = c("Method","Estimate","SE","lowerCI95","upperCI95","pvalue")
  res3[,exposure := "BAT"]
  res3[,outcome := "pineapple"]
  res3[,setting_sex := "men"]
  res3[,setting_Instruments := mySetting]
  res3[,NR_SNPs_total := dim(data_X)[1]]
  res3
}

MR_otherMethods = rbindlist(dumTab3)

#' # Save results ####
#' ***

WriteXLS(x = c("MRTab_FDR","MR_otherMethods"), 
         ExcelFileName=paste0(path_MR_results,MR_prefix,"/04_MR_steroidHormones_Smelling.xlsx"), 
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
