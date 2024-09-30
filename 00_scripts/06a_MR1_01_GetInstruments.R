#' ---
#' title: "Mendelian Randomization - extract instruments for SHs"
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
#' This is the first MR script. Here, I will get the SNP - exposure association for the first MR (SH --> OP). 
#' 
#' In the meeting with Markus, Franz and Katrin (14/05/2024), the following decision on the instrument selection was made: 
#' 
#' - use only pathway SNPs!
#' - create one list of instruments by
#'    1) Select good instrument per exposure and sex-setting
#'    2) Create union between those genes with good instruments
#'    3) Do some pruning (position based)
#'    4) Sensitivity check: only strong instruments per sex-setting
#'    
#' # Initialize ####
#' ***
#' pipeline_name: 15_MR1_01_GetInstruments.R

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
#' ## Load gene list from KEGG 
#' Pathway hsa00140 
geneList = data.table(read_excel(paste0(path_KEGGgenes,"2024_GeneList_SteroidHormonePathway_KEGG_hsa00140.xlsx"),sheet = 2))
cytoList = data.table(read_excel(paste0(path_KEGGgenes,"2024_GeneList_SteroidHormonePathway_KEGG_hsa00140.xlsx"),sheet = 3))

#' remove the gene at chromosome X (sniffin stick data only available on chr 1-22)
cytoList[chr==23,]
cytoList = cytoList[chr<23,]

#' ## Load summary statistics
#' Ruth et al. (2020) - downloaded from GWAS Catalog
ToDoList = fread(paste0(path_SHdata,"2020_Ruth_SexHormones_UKBB_overview.txt"))

statistics = list.files(path = path_SHdata, pattern = ".h.tsv.gz")
statistics2 = unlist(strsplit(statistics,"-"))
statistics2 = statistics2[grepl("GCST",statistics2)]
matched = match(ToDoList$StudyAccession, statistics2)
table(ToDoList$StudyAccession == statistics2[matched])
ToDoList[,statistics := statistics[matched]]
ToDoList[,flag := c(rep("T",3),"E2",rep("SHBG",6),rep("T",3))]

myNames_old = c("hormone","setting",
                "hm_rsid","hm_chrom","hm_pos","hm_effect_allele","hm_other_allele","hm_effect_allele_frequency",
                "hm_beta","standard_error","p_value","zscore","sampleSize")
myNames_new = c("hormone","setting",
                "rsID","CHR","POS","EA","OA","EAF",
                "beta","SE","pvalue","zscore","sampleSize")

dumTab1 = foreach(i=1:dim(ToDoList)[1])%do%{
  #i=1
  myTime = Sys.time()
  myRow = ToDoList[i,]
  message("\nWorking on phenotype ",myRow$Trait," in ",myRow$Setting," (",i," of ",dim(ToDoList)[1],")")
  
  # load data
  myfn1 = paste0(path_SHdata,myRow$statistics)
  erg1 = fread(myfn1)
  erg1[,hormone := myRow$Trait]
  erg1[,setting := myRow$Setting]
  erg1[,sampleSize := myRow$SampleSize]
  erg1[,zscore := hm_beta/standard_error]
  
  stopifnot(myNames_old %in% names(erg1))
  colsOut<-setdiff(colnames(erg1),myNames_old)
  erg1[,get("colsOut"):=NULL]
  setcolorder(erg1,myNames_old)
  names(erg1) = myNames_new
  
  # get list 3: all SNPs of the pathway gene regions
  dumTab2 = foreach(j = 1:dim(cytoList)[1])%do%{
    #j=1
    message("Working on locus number ",j)
    erg4 = copy(erg1)
    myLocus = copy(cytoList)
    myLocus = myLocus[j,]
    
    erg4 = erg4[CHR == myLocus$chr,]
    erg4 = erg4[POS >= myLocus$gene_start-250000,]
    erg4 = erg4[POS <= myLocus$gene_end+250000,]
    
    erg4[,cyto:= myLocus$cytoband]
    erg4
  }
  erg3 = rbindlist(dumTab2)
  erg3[,flag := "PathwaySNPs_all"]

  erg3
  
}
myTab_complete = rbindlist(dumTab1)
save(myTab_complete,file = paste0(path_MR_results,MR_prefix,"/01_SH_instruments_unfiltered.RData"))

myTab_complete[hormone == "BioTestosterone",hormone := "BAT"]
myTab_complete[hormone == "Testosterone",hormone := "TT"]
myTab_complete[hormone == "Estradiol",hormone := "E2"]
myTab_complete[,MAF := EAF]
myTab_complete[EAF>0.5,MAF := 1-EAF]

#' # Get loci per hormone ####
#' ***
myHormones = unique(myTab_complete$hormone)
myHormones = myHormones[-2]

dumTab2 = foreach(i=1:length(myHormones))%do%{
  #i=1
  message("Working on hormone: ",myHormones[i])
  
  myTab = copy(myTab_complete)
  myTab = myTab[hormone == myHormones[i],]
  
  cytoList2 = myTab[abs(zscore)>5,.N,by = c("cyto","setting")]
  cytoList3 = dcast(cytoList2,cyto ~ setting,value.var = "N")
  setDT(cytoList3)
  
  dumTab3 = foreach(j=1:dim(cytoList3)[1])%do%{
    #j=1
    myTab_M = copy(myTab)[setting == "Men" & cyto == cytoList3[j,cyto],]
    myTab_W = copy(myTab)[setting == "Women" & cyto == cytoList3[j,cyto],]
    
    dups = myTab_W[duplicated(rsID),rsID] 
    myTab_W = myTab_W[!is.element(rsID,dups)]
    dups = myTab_M[duplicated(rsID),rsID] 
    myTab_M = myTab_M[!is.element(rsID,dups)]
    
    myTab_M = myTab_M[rsID %in% myTab_W$rsID]
    myTab_W = myTab_W[rsID %in% myTab_M$rsID]
    
    dataset_women = list(beta=myTab_W[,beta], 
                         varbeta=(myTab_W[,SE])^2,
                         N=myTab_W[,sampleSize],
                         snp=myTab_W[,rsID],
                         type="quant",
                         MAF = myTab_W[,MAF],
                         position = myTab_W[,POS])
    dataset_men = list(beta=myTab_M[,beta], 
                         varbeta=(myTab_M[,SE])^2,
                         N=myTab_M[,sampleSize],
                         snp=myTab_M[,rsID],
                         type="quant",
                         MAF = myTab_M[,MAF],
                         position = myTab_M[,POS])
    
    my_res1 = coloc::coloc.abf(dataset1=dataset_women,dataset2=dataset_men)
    my_res2 = my_res1$summary
    trait1 = paste0(myHormones[i]," - Women")
    trait2 = paste0(myHormones[i]," - Men")
    
    # plot data
    par(mar=c(5,4,4,4))
    myXlab=paste0("chromosome ",unique(myTab_M$CHR), " (bp in Mb)")
    myTab_W[,POS2 := POS / 1000000]
    myTab_M[,POS2 := POS / 1000000]

    myTab_W[,plot(POS2, abs(zscore),col = rgb(0,0,1,0.3),pch=19, xlab = myXlab, ylab = "")]
    axis(side = 2, col = "blue", col.ticks = "blue",col.axis="blue")
    mtext(side = 2, line = 2.1, bquote(abs(italic(z-score))), col = "blue", cex = 0.7)

    par(new = T)
    myTab_M[,plot(POS, abs(zscore),axes=F, xlab=NA, ylab=NA, col = rgb(1,0,0,0.2), pch = 17)]
    axis(side = 4, col = "red", col.ticks = "red",col.axis="red")
    mtext(side = 4, line = 2.1, bquote(abs(italic(z-score))), col = "red", cex = 0.7)

    title(paste0(trait1," vs. ",trait2," at ", cytoList3[j,cyto]), cex.main = 1)

    mtext(paste0("H0 | H1 | ... | H4:  ", paste(formatC(round(my_res2[2:6],3), digits=3, format="f" ), collapse = " | ")),cex = 0.7)
    
    x2<-as.data.table(my_res2)
    x3<-t(x2)
    x4<-as.data.table(x3)
    names(x4)<-names(my_res2)
    x4[,cytoband:=cytoList3[j,cyto]]
    x4[,trait1:=trait1]
    x4[,trait2:=trait2]
    snp1 = myTab_W[abs(zscore) == max(abs(zscore)),rsID]
    x4[,bestSNP_trait1 := snp1[1]]
    x4[,bestSNP_trait1_zscore := myTab_W[rsID == snp1[1],abs(zscore)]]
    snp2 = myTab_M[abs(zscore) == max(abs(zscore)),rsID]
    x4[,bestSNP_trait2 := snp2[1]]
    x4[,bestSNP_trait2_zscore := myTab_M[rsID == snp2[1],abs(zscore)]]
    x4

  }
  coloc_res = rbindlist(dumTab3)
  stopifnot(coloc_res$cytoband == cytoList3$cyto)
  coloc_res[PP.H4.abf>0.5,comment := "shared causal SNP"]
  coloc_res[PP.H3.abf>0.5,comment := "two independent SNPs"]
  coloc_res[PP.H2.abf>0.5,comment := "only in men"]
  coloc_res[PP.H1.abf>0.5,comment := "only in women"]
  table(coloc_res$comment)
  coloc_res = cbind(coloc_res,cytoList3[,c(2:4),with=F])
  coloc_res
}
myColoc = rbindlist(dumTab2) 

dumID1 = myColoc[PP.H1.abf>0.5,paste(gsub(" - .*","",trait1),bestSNP_trait1,sep = "_")]
dumID2 = myColoc[PP.H2.abf>0.5,paste(gsub(" - .*","",trait2),bestSNP_trait2,sep = "_")]
dumID31 = myColoc[PP.H3.abf>0.5,paste(gsub(" - .*","",trait1),bestSNP_trait1,sep = "_")]
dumID32 = myColoc[PP.H3.abf>0.5,paste(gsub(" - .*","",trait2),bestSNP_trait2,sep = "_")]
dumID41 = myColoc[PP.H4.abf>0.5 & bestSNP_trait2_zscore<bestSNP_trait1_zscore,
                  paste(gsub(" - .*","",trait1),bestSNP_trait1,sep = "_")]
dumID42 = myColoc[PP.H4.abf>0.5 & bestSNP_trait2_zscore>bestSNP_trait1_zscore,
                  paste(gsub(" - .*","",trait2),bestSNP_trait2,sep = "_")]

dumIDs = c(dumID1,dumID2,dumID31,dumID32,dumID41,dumID42)

#' Coloc was performed for TT, BAT and SHBG, because data from men and women are available. E2 is only availabe in men, hence I can simply pick the best SNP per cytoband
#' 
myTab_E2 = copy(myTab_complete)
myTab_E2 = myTab_E2[hormone == "E2",]
myTab_E2[,absZ := abs(zscore)]
setorder(myTab_E2,-absZ)
myTab_E2 = myTab_E2[!duplicated(cyto),]
myTab_E2 = myTab_E2[absZ>5,]
myTab_E2[,dumID := paste(hormone,rsID,sep="_")]

dumIDs = c(dumIDs,myTab_E2$dumID)

#' Now I can filter for those selected instruments
myTab_complete[,dumID := paste(hormone,rsID,sep="_")]

myTab_filtered = copy(myTab_complete)
myTab_filtered = myTab_filtered[dumID %in% dumIDs,]
myTab_filtered[dumID %in% dumID1, flag := "only in women"]
myTab_filtered[dumID %in% dumID2, flag := "only in men"]
myTab_filtered[dumID %in% dumID31 | dumID %in% dumID32, flag := "two independent SNPs"]
myTab_filtered[dumID %in% dumID41 | dumID %in% dumID42, flag := "shared causal SNP"]
myTab_filtered[hormone == "E2",flag := "only data for men"]

save(myTab_filtered,file = paste0(path_MR_results,MR_prefix,"/01_SNP_on_steroidHormones.RData"))
save(myColoc,file = paste0(path_MR_results,MR_prefix,"/01_Coloc_Women_vs_Men.RData"))

#' # SNP List - updated ####
#' ***
#' In the next script, I have to load and filter all 12 x 3 OP data sets. I need a SNPList to do so. 
#' 
myTab = copy(myTab_filtered)
myTab = myTab[!duplicated(rsID),]
myTab = myTab[,c(3:8)]
myTab[,CHR := as.numeric(CHR)]
setorder(myTab,CHR,POS)

save(myTab,file = paste0(path_MR_results,MR_prefix,"/01_SNPList_instruments.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
