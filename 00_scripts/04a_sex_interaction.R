#' ---
#' title: 'Sex Interaction'
#' subtitle: ''
#' author: 'Andreas, Janne Pott'
#' date: 'Last compiled on `r format(Sys.time(), '%d %B, %Y')`'
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
#' I want to check for sex interaction.
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory

# Plotting of beta-beta plots for each phenotype and one plot for all phenotypes
# combined
# pipeline_name: 11a_sex_interaction.R

#####################################################################################################
##### SEX-SPECIFIC ASSOCIATION ANALYSIS FOR CANDIDATE SNPS FROM CKDGEN GWAS META ANALYSIS: EGFR #####
#####################################################################################################

##### Setup
#

rm(list = ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
library(ggmagnify) # for cropping out a region of the plot
setwd(projectpath)

# create a list of the phenotypes to analyse (only smell without sex subgroup)
phenotype_list = list.files(path_data, ".gz")
phenotype_list = unique(str_match(phenotype_list, "GWASMA_([^_]+)_[^_]+_\\d")[,2])

# phenotype_list = phenotype_list[length(phenotype_list)] #for testing
p_threshold = 5*10^-8 #TODO set threshold (was origninally 5*10^-8)
locus_definition_file = "locus_definition_rsID.csv" #TODO check locus definition file (needs rsID for gwsig)

##### Common Base Data
locus = fread(
  paste0(path_locus_definition, locus_definition_file), 
  dec = ","
)


locus = locus[, c("region", "markerID", "phenotype", "chrom", "pos")]
setnames(
  locus, c("markerID", "chrom", "pos"),
  c("rsID", "chromosome", "position"),
  skip_absent = T
)

setorder(locus, region, chromosome, position)

##### Analysis
j = 0
for (phenotype in phenotype_list) {

  meta_file_all = list.files(
    paste0(path_data),
    pattern = paste0("GWASMA_", phenotype, "_all.*\\.gz")
  )
  meta_file_male = list.files(
    paste0(path_data),
    pattern = paste0("GWASMA_", phenotype, "_male.*\\.gz")
  )
  meta_file_female = list.files(
    paste0(path_data),
    pattern = paste0("GWASMA_", phenotype, "_female.*\\.gz")
  )

  meta_all = fread(
    paste0(path_data, meta_file_all),
    nThread = 30
  )
  meta_male = fread(
    paste0(path_data, meta_file_male),
    nThread = 30
  )
  meta_female = fread(
    paste0(path_data, meta_file_female),
    nThread = 30
  )

  data_list = list(meta_all, meta_male, meta_female)
  lapply(
    data_list, function(x) setnames(
      x,
      c("pFEM", "markerID", "betaFEM", "seFEM"),
      c("P", "rsID", "beta", "SE"),
      skip_absent = T
    )
  )

  ### Unite Results of Males and Females
  meta_all_gwsig = meta_all[(meta_all$P <= p_threshold), ]
  meta_male_gwsig = meta_male[(meta_male$P <= p_threshold), ]
  meta_female_gwsig = meta_female[(meta_female$P <= p_threshold), ]

  union = c(meta_all_gwsig$rsID, meta_male_gwsig$rsID, meta_female_gwsig$rsID)
  union_unique = unique(union)

  meta_male_union = meta_male[match(union_unique, meta_male$rsID),]
  meta_female_union = meta_female[match(union_unique, meta_female$rsID),]

  ### Create Table
  
  # set phenotype to empty vector to prevent errors in dataframe creation if
  # no genomewide significant locus can be found
  if(length(union_unique)==nrow(meta_male_union) && length(union_unique)==nrow(meta_female_union) && length(union_unique)==0){
    phenotype_df = character()
  } else {
    phenotype_df = phenotype
  }
  
  df = data.frame(
    union_unique, phenotype_df, meta_male_union$beta, meta_male_union$SE, meta_male_union$P,
    meta_female_union$beta, meta_female_union$SE, meta_female_union$P, stringsAsFactors = FALSE
  )

  df$mean_diff = rep(NA, nrow(df))
  df$se_diff = rep(NA, nrow(df))
  df$p_diff = rep(NA, nrow(df))

  names(df) = c(
    "SNP", "Trait", "Male_Beta", "Male_SE", "Male_P", "Female_Beta", "Female_SE",
    "Female_P", "mean_diff", "se_diff", "p_diff"
  )

  union_r = c(meta_male$rsID, meta_female$rsID)
  union_unique_r = unique(union_r)
  r = cor(
    meta_male$beta[match(union_unique_r, meta_male$rsID)],
    meta_female$beta[match(union_unique_r, meta_female$rsID)],
    use = "complete.obs", method = "spearman"
  )
  for (i in 1:nrow(df)) {

    z_i = (df$Male_Beta[i] - df$Female_Beta[i])/sqrt(df$Male_SE[i]^2 + df$Female_SE[i]^2 - 2 * r * df$Male_SE[i] * df$Female_SE[i])
    p_i = 2 * pnorm(
      abs(z_i),
      lower.tail = FALSE
    )
    diff_i = (df$Male_Beta[i] - df$Female_Beta[i])
    diff_se_i = sqrt(df$Male_SE[i]^2 + df$Female_SE[i]^2 - 2 * r * df$Male_SE[i] * df$Female_SE[i])

    df$mean_diff[i] = diff_i
    df$se_diff[i] = diff_se_i
    df$p_diff[i] = p_i

  }

  # Match Top-SNPs
  df = df[na.omit(match(locus$rsID, df$SNP)),]

  df$fdr_diff = p.adjust(df$p_diff, method = "fdr")

  col = rep("black", nrow(df))
  col[(df$p_diff <= 0.05)] = "red"


  ### Save data
  write.table(
    df, file = paste0("11_sex_interaction/sex_ia_", phenotype, ".txt"),
    quote = FALSE, sep = "\t", row.names = FALSE
  )

  j = j + 1
  message(str_glue("Done file {j} of {length(phenotype_list)}...\n"))

  

}

# Version for publication -------------------------------------------------
# rsID as label
# only SNPs with sig. difference are labeled
# SNPs that are not significant after FDR correction are not labelelled
# only data for top-phenotype is shown

if(p_threshold == 5e-8){
  total_plot_data = foreach(phenotype = phenotype_list) %do% {
    plot_data = fread(paste0("11_sex_interaction/sex_ia_", phenotype, ".txt"))
    plot_data
  }
  
  sexIA = rbindlist(total_plot_data, fill = T)
  
  #add the region id
  sexIA = merge(sexIA, locus[,c("region", "rsID")], by.x = "SNP", by.y = "rsID")
  
  #TODO row with test for score in region 8 is removed so only the top-phenotype is tested
  # region 7 is removed as it is not independent and therefore should not be used in 
  # multiple testing
  sexIA = sexIA[-c(1,4)]#TODO check
  
  setDT(sexIA)
  sexIA[,fdr_diff := p.adjust(p_diff, method = 'fdr')]
  
  myPlotData<-data.table(rs_id=sexIA$SNP,
                         trait=sexIA$Trait,
                         beta_male = sexIA$Male_Beta,
                         beta_female = sexIA$Female_Beta,
                         se_male = sexIA$Male_SE,
                         se_female = sexIA$Female_SE,
                         meandiff = sexIA$mean_diff,
                         se_meandiff = sexIA$se_diff,
                         meandiff_p = sexIA$p_diff,
                         meandiff_p_FDR = sexIA$fdr_diff,
                         region = sexIA$region)
  myPlotData[meandiff_p>=0.05,sig:='no']
  myPlotData[meandiff_p<0.05,sig:='yes'] # SNPs that are not sig. after FDR correction are considered not sig
  myPlotData[meandiff_p_FDR<0.05,sig:='yes (FDR 5%)']
  myPlotData[,gene2:='']
  # myPlotData[meandiff_p_FDR<0.05,gene2:=gsub(':.*','',rs_id)]
  myPlotData[,gene2:=paste0(region)]
  myPlotData[,gene3:='']
  myPlotData[meandiff_p<0.05,gene3:=gsub(':.*','',rs_id)]
  myPlotData[, sex.higherEffect := 'none']
  myPlotData[(abs(beta_female) > abs(beta_male)) & sig == 'yes (FDR 5%)', sex.higherEffect := 'female']
  myPlotData[(abs(beta_male) > abs(beta_female)) & sig == 'yes (FDR 5%)', sex.higherEffect := 'male']
  
  #set colours
  cols = c("none" = '#000000', "female" = '#B2182B', "male" = '#2166AC')
  
  
  
  
  
  # create rsID as label for SNPs with sig. interaction
  locus_def = fread(paste0(path_locus_definition, "locus_definition_rsID.csv"),
                    dec = ",") #TODO check locus definition file with rsIDs
  locus_def = locus_def[,c("markerID", "rsID")]
  myPlotData = merge(myPlotData, locus_def, by.x = "rs_id", by.y = "markerID", all = F)
  
  
  #select snps for manual labeling
  snp1_label = subset(myPlotData, rsID == "rs116058752")
  snp2_label = subset(myPlotData, rsID == "rs56320200")
  snp3_label = subset(myPlotData, rsID == "rs61902559")
  
  myPlot1 = ggplot(myPlotData,
                   aes(x = beta_male, y = beta_female, color = sex.higherEffect)) +
    # facet_wrap(~trait, scales = 'free') +
    ylim(c(-2, 1.5)) + #aa included
    # ylim(c(-2.1, 1.5)) +
    xlim(c(-1.5, 0.75)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = 'grey',
      linetype = 'dashed',
      size = 1.25
    ) +
    geom_hline(
      yintercept = 0,
      color = 'grey',
      linetype = 'dashed',
      size = 1.15
    ) +
    geom_vline(
      xintercept = 0,
      color = 'grey',
      linetype = 'dashed',
      size = 1.15
    ) +
    geom_point(size = 3) +
    geom_errorbar(aes(
      ymin = beta_female - 1.96 * se_female,
      ymax = beta_female + 1.96 * se_female
    ),
    width = 0.02) +
    geom_errorbarh(aes(xmin = beta_male - 1.96 * se_male, xmax = beta_male + 1.96 *
                         se_male),
                   height = 0.02) +
    theme_bw(base_size = 10) +
    scale_colour_manual(values = cols) +
    theme(
      plot.title = element_text(hjust = 0, size = 22, face = 'bold'),
      axis.title.x = element_text(size = 22, face = 'bold'),
      axis.title.y = element_text(size = 22, face = 'bold'),
      axis.text = element_text(size = 20, face = 'bold'),
      strip.text = element_text(size = 20)
    ) +
    labs(x = 'Effect Size Male',
         y = 'Effect Size Female', color = 'SNPs with \\nsex interaction') +
    geom_label_repel(
      data = snp1_label,
      aes(x = beta_male, y = beta_female, label = rsID),
      nudge_x = 0.3,
      nudge_y = 0.1,
      max.overlaps = Inf,
      size = 7
    ) +
    geom_label_repel(
      data = snp2_label,
      aes(x = beta_male, y = beta_female, label = rsID),
      nudge_x = 0.3,
      nudge_y = 0.1,
      max.overlaps = Inf,
      size = 7
    ) +
    geom_label_repel(
      data = snp3_label,
      aes(x = beta_male, y = beta_female, label = rsID),
      nudge_x = -0.3,
      nudge_y = 0.1,
      max.overlaps = Inf,
      size = 7
    ) +
    guides(label = 'none', color = 'none') 
    # geom_magnify(
    #   # from = c(-0.01, 0.6, -0.2, 0.55),
    #   #values for aa_included
    #   from = c(-0.43, 0.15, -0.51 ,0.07), # aa excluded
    #   #crop out dense regions
    #   # to = c(0, 1.5, -2, -0.52),
    #   #values for aa_included
    #   to = c(-2.7, -0.5, 0.22, 2), # aa excluded
    #   shadow = F,
    #   recompute = T,
    #   # aspect = "fixed",
    #   # proj = "single"
    # )
  
  # print(myPlot1)
  
  png(filename = paste0('11_sex_interaction/BetaBeta_sexIA_publication.png'),
      width = 6000, height = 6000, res=600)
  print(myPlot1)
  dev.off()
  
  #save new test data because only the top-phenotype is tested and therefore the 
  # FRD corrected p-values change
  setnames(myPlotData, "gene2", "region", skip_absent = T)
  setkey(myPlotData, region)
  write.table(
    myPlotData,
    "11_sex_interaction/sex_ia.csv",
    col.names = T,
    row.names = F,
    quote = F,
    sep = ";",
    dec = ","
  )
}
# Unlabelled version for suggestive hits -----------------------------------

if(p_threshold > 5e-8){
  total_plot_data = foreach(phenotype = phenotype_list) %do% {
    plot_data = fread(paste0("11_sex_interaction/sex_ia_", phenotype, ".txt"))
    plot_data
  }
  
  sexIA = rbindlist(total_plot_data, fill = T)
  
  #add the region id
  sexIA = merge(sexIA, locus[,c("region", "rsID")], by.x = "SNP", by.y = "rsID")
  
  
  setDT(sexIA)
  sexIA[,fdr_diff := p.adjust(p_diff, method = 'fdr')]
  
  myPlotData<-data.table(rs_id=sexIA$SNP,
                         trait=sexIA$Trait,
                         beta_male = sexIA$Male_Beta,
                         beta_female = sexIA$Female_Beta,
                         se_male = sexIA$Male_SE,
                         se_female = sexIA$Female_SE,
                         meandiff = sexIA$mean_diff,
                         se_meandiff = sexIA$se_diff,
                         meandiff_p = sexIA$p_diff,
                         meandiff_p_FDR = sexIA$fdr_diff,
                         region = sexIA$region)
  myPlotData[meandiff_p>=0.05,sig:='no']
  myPlotData[meandiff_p<0.05,sig:='yes']
  myPlotData[meandiff_p_FDR<0.05,sig:='yes (FDR 5%)']
  myPlotData[,gene2:='']
  # myPlotData[meandiff_p_FDR<0.05,gene2:=gsub(':.*','',rs_id)]
  myPlotData[,gene2:=paste0(region)]
  myPlotData[,gene3:='']
  myPlotData[meandiff_p<0.05,gene3:=gsub(':.*','',rs_id)]
  myPlotData[, sex.higherEffect := 'none']
  myPlotData[(abs(beta_female) > abs(beta_male)) & sig == 'yes (FDR 5%)', sex.higherEffect := 'female']
  myPlotData[(abs(beta_male) > abs(beta_female)) & sig == 'yes (FDR 5%)', sex.higherEffect := 'male']
  
  #set colours
  cols = c("none" = '#000000', "female" = '#B2182B', "male" = '#2166AC')
  
  left_label = subset(myPlotData, beta_male<0)
  right_label = subset(myPlotData, beta_male>0)
  
  myPlot1 = ggplot(myPlotData,
                   aes(x = beta_male, y = beta_female, color = sex.higherEffect)) +
    # facet_wrap(~trait, scales = 'free') +
    ylim(c(-2, 2)) + #aa included
    xlim(c(-2, 2)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = 'grey',
      linetype = 'dashed',
      size = 1.25
    ) +
    geom_hline(
      yintercept = 0,
      color = 'grey',
      linetype = 'dashed',
      size = 1.15
    ) +
    geom_vline(
      xintercept = 0,
      color = 'grey',
      linetype = 'dashed',
      size = 1.15
    ) +
    geom_point(size = 3) +
    geom_errorbar(aes(
      ymin = beta_female - 1.96 * se_female,
      ymax = beta_female + 1.96 * se_female
    ),
    width = 0.02) +
    geom_errorbarh(aes(xmin = beta_male - 1.96 * se_male, xmax = beta_male + 1.96 *
                         se_male),
                   height = 0.02) +
    theme_bw(base_size = 10) +
    scale_colour_manual(values = cols) +
    theme(
      plot.title = element_text(hjust = 0, size = 22, face = 'bold'),
      axis.title.x = element_text(size = 22, face = 'bold'),
      axis.title.y = element_text(size = 22, face = 'bold'),
      axis.text = element_text(size = 20, face = 'bold'),
      strip.text = element_text(size = 20)
    ) +
    labs(x = 'Effect Size Male',
         y = 'Effect Size Female', color = 'SNPs with \\nsex interaction') +
    guides(label = 'none', color = 'none')
  
  # print(myPlot1)
  
  png(filename = paste0('11_sex_interaction/BetaBeta_sexIA_suggestive.png'),
      width = 6000, height = 6000, res=600)
  print(myPlot1)
  dev.off()
  
  setnames(myPlotData, "gene2", "order_region", skip_absent = T)
  myPlotData[,order_region := as.numeric(str_remove(order_region, "S"))]
  setkey(myPlotData, order_region)
  
  write.table(
    myPlotData,
    "11_sex_interaction/sex_ia_suggestive.csv",
    col.names = T,
    row.names = F,
    quote = F,
    sep = ";",
    dec = ","
  )
  #' # Session Info ####
  #' ***
}

sessionInfo()
message(
  "\nTOTAL TIME : ", round(
    difftime(Sys.time(), time0, units = "mins"),
    3
  ),
  " minutes"
)


