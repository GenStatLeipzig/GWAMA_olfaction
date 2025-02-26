#' ---
#' title: "Mendelian Randomization summary"
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
#' This is the last script for the MR analyses. Here, I just want to generate one excel file with all relevant tables.
#'
#' pipeline_name: 15_MR_summary.R
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(paste0(projectpath, "00_scripts/"))
.libPaths()

MR1_prefix = "MR1_steroidHormone_on_Smelling"
MR2_prefix = "MR2_Smelling_on_neuroDiseases"
# MR3_prefix = "" # deprecated analysis
MR4_prefix = "MR3_neuroDiseases_on_Smelling"

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(
    Table = paste0("S", c(1:7)),
    Title = c(
      "MR1 - hormone --> stick; summary statistics (data used)",
      "MR1 - hormone --> stick; MR-IVW results (hier. FDR corrected)",
      "MR2 - stick --> disease; summary statistics (data used)",
      "MR2 - stick --> disease; Wald ratios or MR-IVW results (hier. FDR corrected)",
      "MR4 - disease --> stick; summary statistics (data used)",
      "MR4 - disease --> stick; MR-IVW results (hier. FDR corrected)",
      "MR4 - disease --> stick; other MR-models for significant combination"
    )
  )

  tab0
}

#' # Get table S1 ####
#' ***
#' What do I want in here?
#'
#' - SNP: cyto, gene, rsID, chr, pos, EA, OA
#' - exposure: hormone, EAF, beta, se, pval, zscore, sample size
#' - outcome: stick, EAF, beta, se, pval, zscore, sample size
#'
{
  load(paste0(path_MR_results, MR1_prefix, "/01_SNP_on_steroidHormones.RData"))
  load(paste0(path_MR_results, MR1_prefix, "/02_SNP_on_olfactoryPerception.RData"))

  test = myTab_filtered[, max(abs(zscore)), by = rsID]
  test = test[V1 < 5.45, ]
  myTab_filtered = myTab_filtered[!is.element(rsID, test$rsID)]

  cytoList = data.table(read_excel(paste0(path_KEGGgenes, "2024_GeneList_SteroidHormonePathway_KEGG_hsa00140.xlsx"), sheet = 3))
  table(is.element(myTab_filtered$cyto, cytoList$cytoband))
  matched = match(myTab_filtered$cyto, cytoList$cytoband)
  myTab_filtered[, genes := cytoList[matched, genes]]

  myTab_filtered[, exposure := paste(hormone, setting, sep = "_")]
  myExposures = unique(myTab_filtered$exposure)
  myTab_OP[, outcome := paste(phenotype, setting, sep = "_")]
  myOutcomes = unique(myTab_OP$outcome)

  dumTab1 = foreach(i = 1:length(myExposures)) %do% {
    # i=1
    data_X = copy(myTab_filtered)
    data_X = data_X[exposure == myExposures[i], ]

    dumTab2 = foreach(j = 1:length(myOutcomes)) %do% {
      # j=1
      data_Y = copy(myTab_OP)
      data_Y = data_Y[outcome == myOutcomes[j]]

      # match data
      matched = match(data_X$rsID, data_Y$rsID)
      data_Y = data_Y[matched, ]
      filt = is.na(data_Y$beta) | is.na(data_X$beta)
      data_Y = data_Y[!filt, ]
      data_X = data_X[!filt, ]
      stopifnot(data_X$rsID == data_Y$rsID)

      # get nice output
      data = cbind(data_X[, c(14, 18, 3:7, 19, 8:13)], data_Y[, c(19, 8, 12:16)], data_X[, flag])
      names(data)[9:14] = paste0("exposure_", names(data)[9:14])
      names(data)[16:21] = paste0("outcome_", names(data)[16:21])
      data
    }
    data1 = rbindlist(dumTab2)
    data1
  }
  tab1 = rbindlist(dumTab1)
  tab1 = tab1[grepl("Comb", exposure) & grepl("all", outcome) |
    grepl("_Women", exposure) & grepl("female", outcome) |
    grepl("_Men", exposure) & grepl("_male", outcome)]
  tab1[, sex := "all"]
  tab1[grepl("_Women", exposure), sex := "women"]
  tab1[grepl("_Men", exposure), sex := "men"]
  tab1[, exposure := gsub("_.*", "", exposure)]
  tab1[, outcome := gsub("_.*", "", outcome)]
  tab1 = tab1[, c(1:7, 23, 8:22)]

  # add flag to indicate in which setting each SNP was used (main: strong SNPs, sens1: same SNPs, sens2: colocalizing SNPs)
  tab1[, main := F]
  tab1[abs(exposure_zscore) > 5.45, main := T]
  tab1[, sens1 := T]
  tab1[, sens2 := F]
  tab1[grepl("shared", V3), sens2 := T]
  tab1[, table(sens2, main)]
  tab1[, V3 := NULL]
  tab1
}

#' # Get table S2-S3 ####
#' ***
#' MR1: hormone levels on smell recognition
{
  loaded = load(paste0(path_MR_results, MR1_prefix, "/04_MRresults_IVW_FDRadj.RData"))
  tab2 = get(loaded)
  names(tab2)
  tab2[, method := "IVW"]
  tab2[, sex := "all"]
  tab2[grepl("Women", exposure), sex := "women"]
  tab2[grepl("Men", exposure), sex := "men"]
  tab2 = tab2[, c(13, 1, 2, 14, 10, 3, 4:7, 11, 12, 8, 9)]
  tab2[, exposure := gsub("_BMIadj", "BMIadj", exposure)]
  tab2[, exposure := gsub("_.*", "", exposure)]
  setnames(tab2, "SNPset", "mode")
  setnames(tab2, "NR_SNPs_total", "SNPs")
  tab2[mode == "strongSNPs", mode := "main"]
  tab2[mode == "sameSNPs", mode := "sens1"]
  tab2[mode == "sameStrongSNPs", mode := "sens2"]

  setnames(tab2, "beta_IVW", "Estimate")
  setnames(tab2, "SE_IVW", "SE")
  setnames(tab2, "pval_IVW", "pvalue")
  setnames(tab2, "pval_adjusted", "pvalue_adj1")
  setorder(tab2, mode)

  tab2[exposure == "E2"]
  tab2[exposure == "E2", mode := "main"]

  myModes = unique(tab2$mode)
  mySexes = unique(tab2$sex)

  dumTabMR1_1 = foreach(i = 1:length(myModes)) %do% {
    # i=1
    dumTabMR1_2 = foreach(j = 1:length(mySexes)) %do% {
      # j=2
      data = copy(tab2)
      data = data[mode == myModes[i] & sex == mySexes[j]]

      data[, FDR1 := p.adjust(p = pvalue, method = "fdr"), by = outcome]
      data[, FDR2 := p.adjust(p = FDR1, method = "fdr"), by = exposure]
      data[, hierFDR := FDR2 < 0.05]

      data = data[, c(1:9, 15, 16, 12, 10, 13, 14)]
      data
    }
    dumTabMR1_2 = rbindlist(dumTabMR1_2)
    dumTabMR1_2
  }
  tab2 = rbindlist(dumTabMR1_1)

  # tab3 can be ignored, as nothing is significant - no reason to run the other methods!
  loaded = load(paste0(path_MR_results, MR1_prefix, "/04_MRresults_otherMethods.RData"))
  tab3 = get(loaded)
  names(tab3)
  tab3 = tab3[, c(1, 7, 8, 9, 10, 11, 2:6)]
  setnames(tab3, "Method", "method")
  setnames(tab3, "setting_sex", "sex")
  setnames(tab3, "setting_Instruments", "mode")
  setnames(tab3, "NR_SNPs_total", "SNPs")
  tab3[mode == "strongSNPs", mode := "main"]
  tab3[mode == "sameSNPs", mode := "sens1"]
  tab3[mode == "sameStrongSNPs", mode := "sens2"]
  setorder(tab3, mode)
}

#' # Get table S4-S5 ####
#' ***
#' MR2: smell recognition on neuro-degenerative diseases
{
  loaded = load(paste0(path_MR_results, MR2_prefix, "/02_MRbySNP_Smelling_NeuroDiseases.RData"))
  tab4 = get(loaded)
  names(tab4)
  tab4[, method := "WaldRatio"]
  tab4[, mode := "sensitivity"]
  tab4[bestPhenotype == phenotype, mode := "main"]
  tab4[grepl("female", exposure_setting), exposure_setting := "women"]
  tab4[grepl("male", exposure_setting), exposure_setting := "men"]
  setnames(tab4, "exposure_setting", "sex")
  tab4[outcome == "PARKINSON", outcome_sampleSize := 4681 + 1660 + 407500 + 403249]
  tab4[outcome == "PARKINSON", outcome_cases := 4681 + 1660]
  tab4[outcome == "AD_WIDE", outcome_sampleSize := 15617 + 839 + 396564 + 410833]
  tab4[outcome == "AD_WIDE", outcome_cases := 15617 + 839]
  tab4[outcome == "ALZHEIMER", outcome_sampleSize := 10520 + 831 + 401661 + 419700]
  tab4[outcome == "ALZHEIMER", outcome_cases := 10520 + 831]
  tab4 = tab4[, c(26, 27, 9, 3:8, 10:20, 28, 29, 21:23, 25)]
  setorder(tab4, mode)

  loaded = load(paste0(path_MR_results, MR2_prefix, "/03_MRIVW_Smelling_NeuroDiseases.RData"))
  tab5 = get(loaded)
  names(tab5)
  tab5[, method := "IVW"]
  tab5[, sex := "all"]
  tab5[grepl("female", exposure), sex := "women"]
  tab5[grepl("_male", exposure), sex := "men"]
  tab5 = tab5[, c(12, 1, 2, 13, 10, 3, 4:7, 11, 8, 9)]
  tab5[, exposure := gsub("_.*", "", exposure)]
  setnames(tab5, "comment", "mode")
  setnames(tab5, "NR_SNPs_total", "SNPs")
  tab5[, mode := "main"]
  tab5[sex != "all", mode := "sensitivity"]
  setnames(tab5, "beta_IVW", "Estimate")
  setnames(tab5, "SE_IVW", "SE")
  setnames(tab5, "pval_IVW", "pvalue")
  setnames(tab5, "pval_adj", "pvalue_adj")

  # resort tables: main = Wald ratio in case of single IV (with best associated sex), IVW all in case of multiple IVs
  # sensitivity: other Wald ratios (multiple SNPs, other sexes) and IVW other sexes
  tab4_2 = tab4[, c(1, 9, 15, 3, 2, 4, 23, 24, 25)]
  setnames(tab4_2, "rsID", "SNPs")
  tab4_2[, SNPs := 1]
  setnames(tab4_2, "MR_beta", "Estimate")
  setnames(tab4_2, "MR_SE", "SE")
  setnames(tab4_2, "MR_pvalue", "pvalue")
  tab4_2[, Fstat := tab4$exposure_beta / tab4$exposure_SE]
  tab4_2 = tab4_2[exposure %in% c("lemon", "orange")]
  tab5_2 = rbind(tab5, tab4_2, fill = T)

  myModes = unique(tab5_2$mode)
  mySexes = unique(tab5_2$sex)

  dumTabMR1_1 = foreach(i = 1:length(myModes)) %do% {
    # i=1
    if (myModes[i] == "main") {
      data = copy(tab5_2)
      data = data[mode == myModes[i]]

      data[, FDR1 := p.adjust(p = pvalue, method = "fdr"), by = outcome]
      data[, FDR2 := p.adjust(p = FDR1, method = "fdr"), by = exposure]
      data[, hierFDR := FDR2 < 0.05]

      data = data[, c(1:9, 14, 15, 16, 10, 12, 13)]
    } else {
      dumTabMR1_2 = foreach(j = 1:length(mySexes)) %do% {
        # j=1
        data = copy(tab5_2)
        data = data[mode == myModes[i] & sex == mySexes[j]]

        data[, FDR1 := p.adjust(p = pvalue, method = "fdr"), by = outcome]
        data[, FDR2 := p.adjust(p = FDR1, method = "fdr"), by = exposure]
        data[, hierFDR := FDR2 < 0.05]

        data = data[, c(1:9, 14, 15, 16, 10, 12, 13)]
        data
      }
      data = rbindlist(dumTabMR1_2)
    }
    data
  }
  tab5_3 = rbindlist(dumTabMR1_1, use.names = T)
  tab4 = tab4[, c(4:8, 3, 9:25)]
  names(tab4) = gsub("MR_", "WaldRatio_", names(tab4))

  names(tab5_3)
}


#' # Get table S9 ####
#' ***
#' What do I want in here?
#'
#' - SNP: rsID, chr, pos, EA, OA
#' - exposure: disease, EAF, beta, se, pval, zscore, sample size (cases & controls)
#' - outcome: stick, EAF, beta, se, pval, zscore, sample size
#'
{
  load(paste0(path_MR_results, MR4_prefix, "/01_SNP_on_neuroDiseases.RData"))
  load(paste0(path_MR_results, MR4_prefix, "/02_SNP_on_olfactoryPerception.RData"))

  #' check outcome data: I want all SNPs to be available in all outcomes! There is one variant that was filtered in the sex-stratified data due to MAC problems.
  #'
  test = myTab_OP[, .N, by = rsID]
  test[N != 39, ]
  myTab_filtered[rsID == test[N != 39, rsID]]
  myTab_OP = myTab_OP[!is.element(rsID, test[N != 39, rsID]), ]
  myTab_filtered = myTab_filtered[!is.element(rsID, test[N != 39, rsID]), ]

  setnames(myTab_filtered, "disease", "exposure")
  myExposures = unique(myTab_filtered$exposure)
  myTab_OP[, outcome := paste(phenotype, setting, sep = "_")]
  myOutcomes = unique(myTab_OP$outcome)

  dumTab1 = foreach(i = 1:length(myExposures)) %do% {
    # i=1
    data_X = copy(myTab_filtered)
    data_X = data_X[exposure == myExposures[i], ]

    dumTab2 = foreach(j = 1:length(myOutcomes)) %do% {
      # j=1
      data_Y = copy(myTab_OP)
      data_Y = data_Y[outcome == myOutcomes[j]]

      # match data
      matched = match(data_X$rsID, data_Y$rsID)
      data_Y = data_Y[matched, ]
      filt = is.na(data_Y$beta) | is.na(data_X$beta)
      data_Y = data_Y[!filt, ]
      data_X = data_X[!filt, ]
      stopifnot(data_X$rsID == data_Y$rsID)

      # get nice output
      data = cbind(data_X[, c(3, 6:9, 1, 10:18)], data_Y[, c(1, 2, 8, 12:16)])
      names(data)[7:15] = paste0("exposure_", names(data)[7:15])
      names(data)[18:23] = paste0("outcome_", names(data)[18:23])
      names(data)[16] = "outcome"
      data
    }
    data1 = rbindlist(dumTab2)
    data1
  }
  tab9 = rbindlist(dumTab1)
  names(tab9)
  setnames(tab9, "setting", "sex")
  tab9[grepl("female", sex), sex := "women"]
  tab9[grepl("male", sex), sex := "men"]
  tab9 = tab9[, c(1:5, 17, 6:16, 18:23)]
}
tab9 = tab9[outcome == "SCORE"]
#' # Get table S10-S11 ####
#' ***
#' MR4: disease on smell recognition
{
  loaded = load(paste0(path_MR_results, MR4_prefix, "/04_MRresults_IVW_FDRadj.RData"))
  tab10 = get(loaded)
  tab10 = tab10[outcome == "SCORE"]
  names(tab10)
  tab10[, method := "IVW"]
  setnames(tab10, "setting", "sex")
  tab10[grepl("female", sex), sex := "women"]
  tab10[grepl("male", sex), sex := "men"]
  tab10 = tab10[, c(14, 1, 2, 11, 3, 4:7, 12, 8, 9)]
  setnames(tab10, "NR_SNPs_total", "SNPs")

  setnames(tab10, "beta_IVW", "Estimate")
  setnames(tab10, "SE_IVW", "SE")
  setnames(tab10, "pval_IVW", "pvalue")
  setnames(tab10, "pval_adjusted", "pvalue_adj")

  mySexes = unique(tab10$sex)

  dumTabMR4_2 = foreach(j = 1:length(mySexes)) %do% {
    # j=1
    data = copy(tab10)
    data = data[sex == mySexes[j]]

    data[, FDR1 := p.adjust(p = pvalue, method = "fdr"), by = outcome]
    data[, FDR2 := p.adjust(p = FDR1, method = "fdr"), by = exposure]
    data[, hierFDR := FDR2 < 0.05]

    data = data[, c(1:8, 13, 14, 15, 9, 11, 12)]
    data
  }
  tab10 = rbindlist(dumTabMR4_2)

  loaded = load(paste0(path_MR_results, MR4_prefix, "/04_MRresults_otherMethods.RData"))
  tab11 = get(loaded)
  names(tab11)
  tab11 = tab11[, c(1, 7, 8, 9, 10, 2:6)]
  setnames(tab11, "Method", "method")
  setnames(tab11, "setting_sex", "sex")
  setnames(tab11, "NR_SNPs_total", "SNPs")

  # select replications for significant combinations of outcome, exposure and sex
  sig_combs = tab10[hierFDR == TRUE, .(outcome, exposure, sex)]
  sig_combs[sex == "women", sex := "female"]
  sig_combs[sex == "men", sex := "male"]
  for (r in 1:nrow(sig_combs)) {
    tab11[outcome == sig_combs[r, outcome] & exposure == sig_combs[r, exposure] & sex == sig_combs[r, sex], include := TRUE]
  }
  tab11 = tab11[include == TRUE]
  tab11[, include := NULL]
}

#' # Save tables ###
#' ***


tosave4 = data.table(
  data = c(
    "tab0", "tab1", "tab2", "tab4", "tab5_3",
    "tab9", "tab10", "tab11"
  ),
  SheetNames = c(
    "Content", "TableS1", "TableS2", "TableS3", "TableS4",
    "TableS5", "TableS6", "TableS7"
  )
)
excel_fn = paste0(path_MR_results, "/SupplementalTables_MR_FDRupdate.xlsx")
WriteXLS(tosave4$data,
  ExcelFileName = excel_fn,
  SheetNames = tosave4$SheetNames,
  AutoFilter = T,
  BoldHeaderRow = T,
  FreezeRow = 1
)

save(tab0, tab1, tab2, tab3, tab4, tab5_3, tab9, tab10, tab11,
  file = paste0(path_MR_results, "/SupplementalTables_MR_FDRupdate.RData")
)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "mins"), 3), " minutes")
