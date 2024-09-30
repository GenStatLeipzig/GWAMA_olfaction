#' ---
#' title: "Mendelian Randomization 2 - run MR analyses"
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
#' This is the second script for MR2. Here, I will run the MR analyses per exposure (different sniffin sticks) and each outcome (3 neuro diseases).
#'  
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR2_02_MRbySNP.R

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
load(paste0(path_MR_results,MR_prefix,"/01_SNP_on_olfactoryPerception.RData"))
load(paste0(path_MR_results,MR_prefix,"/01_SNP_on_NeuroOutcomes.RData"))

myTab1 = fread(paste0("../", path_locus_definition, "/locus_definition_rsID.csv"), dec = ",")
myNames = c("rsID","markerID","chrom","pos","aa","ea","nWeightedMAF","nWeightedInfoScore",
            "phenotype","betaFEM","seFEM","pFEM","I2","totalN")
stopifnot(myNames %in% names(myTab1))
colsOut<-setdiff(colnames(myTab1),myNames)
myTab1[,get("colsOut"):=NULL]
setcolorder(myTab1,myNames)

matched = match(myTab2$markerID,myTab1$markerID)
myTab2[,bestPhenotype := myTab1[matched,phenotype]]

#' # SNP-wise MR ####
#' ***
#' 
#' - SNP: bestPhenotype, markerID, rsID (if available), chr, pos, EA, OA
#' - exposure: stick, setting, EAF, beta, se, pval, sample size
#' - outcome: outcome, EAF, beta, se, pval
#' - MR: estimate, SE, pval, adj. pval
#' 
myOutcomes = myTab3[,unique(outcome)]
mySettings = c("all","female","male")

dumTab2 = foreach(i = 1:length(myOutcomes))%do%{
  #i=1
  myTab6 = copy(myTab3)
  myTab6 = myTab6[outcome == myOutcomes[i],]
  
  dumID = myTab1$dumID
  
  dumTab4 = foreach(j=1:3)%do%{
    #j=1
    myTab4 = copy(myTab2)
    myTab4 = myTab4[setting == mySettings[j],]
    
    stopifnot(myTab4$markerID == myTab6$matchingID)
    stopifnot(myTab4$EA == myTab6$EA)
    
    test = MRfunction_jp(betaX = myTab4$beta,seX = myTab4$SE,betaY = myTab6$beta,seY = myTab6$SE)
  
    result = cbind(myTab4[,c(17,1:8,12,13,14,16)],
                   myTab6[,c(1,7,8,9,12:14)],test[,c(1,3,5)])
    result = result[,c(1,4,14,5:8,2,3,9:13,20,18,19,15:17,21:23)]
    names(result)[8] = "exposure"
    names(result)[9:14] = paste0("exposure_",names(result)[9:14])
    names(result)[16:20] = paste0("outcome_",names(result)[16:20])
    names(result)[21:23] = c("MR_beta","MR_SE","MR_pvalue")
    
    setnames(result,"outcome_pval","outcome_pvalue")
    result
  }

  test = rbindlist(dumTab4)
  test
}
myMRresults2 = rbindlist(dumTab2)

myMRresults2[,dummy := gsub("_.*","",bestPhenotype)]
myMRresults2 = myMRresults2[dummy == exposure,]
myMRresults2[,dummy := NULL]

myMRresults2[,phenotype:=paste(exposure,exposure_setting,sep="_")]

#' Filter for genome-wide significant instruments only
myMRresults2 = myMRresults2[exposure_pvalue<5e-8,]
myMRresults2[,table(rsID,exposure_setting)]

#' Perform FDR per SNP
myMRresults2[MR_pvalue<0.05,]
myMRresults2[MR_pvalue<0.05 & phenotype == bestPhenotype,]
myMRresults2[,MR_pval_adj := p.adjust(p = MR_pvalue,method = "fdr"),by=markerID]
myMRresults2[MR_pval_adj<0.05,]

#' # Visualize results ####
#' ***
#' I want to make some Forest Plot for the significant single SNP results (lemon all on ALZHEIMER).
#' 
plotData1 = copy(myMRresults2)
myIDs = myMRresults2[MR_pvalue<0.05,unique(markerID)]
plotData1 = plotData1[markerID %in% myIDs,]
setorder(plotData1,MR_beta)

plotData1[,lowerCI95 := MR_beta-1.96*MR_SE]
plotData1[,upperCI95 := MR_beta+1.96*MR_SE]
xmin = min(c(0,plotData1$lowerCI95, plotData1$upperCI95), na.rm = T)
xmax = max(c(0,plotData1$lowerCI95, plotData1$upperCI95), na.rm = T)
xmin2 = -0.75
xmax2 = 0.5
# if(xmax2>5) xmax2 = 5
# if(xmin2<-5) xmin2 = -5
xrange = xmax2 - xmin2
xmin_margin = xmin2 - 0.7*xrange
xmax_margin = xmax2 + 0.7*xrange

bestX = gsub("_.*","",plotData1$bestPhenotype[1])
bestset = gsub(".*_","",plotData1$bestPhenotype[1])

fn = paste0(path_MR_results,MR_prefix,"/02_ForestPlot_SNP_",unique(plotData1$rsID),".png")
png(filename = fn,width = 2000, height = 1200, res=250)
par(mar=c(5,6,0,4))
par(font=1)
dets = forest(x=plotData1[,MR_beta], 
              sei = plotData1[,MR_SE],
              header = paste0("Lemon recognition on neuredegenerative diseases"),
              xlab = paste(unique(plotData1$rsID)," - best phenotype in GWAS: ",bestX, "(",bestset,")"),
              slab = paste(plotData1[,outcome]," - ",plotData1[,exposure_setting]),
              ylim = c(-0, nrow(plotData1)+3) , 
              cex = 1, 
              xlim = c(xmin_margin, xmax_margin), 
              alim = c(xmin2, xmax2),at = c(-0.75,-0.5,-0.25,0,0.25,0.5))
par(font=4)
text(min(dets$xlim), max(dets$ylim-1.5), "Disease - sex-setting", pos=4)
dev.off()


#' # Save results ####
#' ***
#' I want an Excel sheet - what do I want in here? 
#' 
#' - SNP: rsID, locus NR, cyto, chr, pos, EA, OA
#' - exposure: stick, setting, EAF, beta, se, pval
#' - outcome: outcome, EAF, beta, se, pval
#' - MR: estimate, SE, pval, adj. pval
#' 

WriteXLS(myMRresults2, 
         ExcelFileName=paste0(path_MR_results,MR_prefix,"/02_MRbySNP_Smelling_NeuroDiseases.xlsx"), 
         SheetNames="MRbySNP", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(myMRresults2,file = paste0(path_MR_results,MR_prefix,"/02_MRbySNP_Smelling_NeuroDiseases.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
