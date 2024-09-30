# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-05-06
#
# Script Description: Create Manhattan plot for publication with separate plots for 
# combined, males and females
#
#
# Notes: see https://github.com/juliedwhite/miamiplot/blob/master/vignettes/scratch_miamiplots.Rmd
# for code template
# pipeline_name: 04c_mh_plots_all_gw_hits_publication.R
#

# INIT --------------------------------------------------------------------
rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggmagnify))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(lemon))

set.seed(456734)
# VARIABLES ---------------------------------------------------------------

locusdef_file = "locus_definition_rsID.csv"

image_name = "all_gw_hits_publication.png" #TODO check name of saved image
#TODO check name for saved plot data
# data_name = paste0("plot_data_all_gw_hits_aa_excluded_publication", Sys.Date(), ".RData")
data_name = "plot_data_all_gw_hits_aa_excluded_publication_2024-05-07.RData"


# qc filters
maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

chrom_number = 22

max.cores = 40

max_p = 20 #maximal log p value to plot
filter_max = T #check wheter filtering for max values

#set color scheme for significant traits
col_dir = c(
  "grey30" = "grey30",
  "grey" = "grey",
  "cinnamon" = "brown",
  "fish" = "royalblue",
  "lemon" = "yellowgreen",
  "orange" = "orange",
  "pineapple" = "darkgoldenrod4",
  "score" = "black"
  
)

# additional annotation of loci (list in order of loci)
cand_genes = c("\nGSX2, FIP1L", "\nADCY2", "\nFBXL17", "\nOR cluster", "\nTAAR5", "\nOR cluster", "", "\nOR cluster", "\nOR cluster", "", "\nOR cluster")

sex_ia = rep(F, 11)
sex_ia[c(2,9,10)] = T

novelty = c(T,T,T,T,F,F,F,T,F,T,T) #locus 7 is set to known as we do not want to hightlight the not independent locus
dependent_regions = c(7) #add list of loci that are not considered independent


collect_plot_data = T #create new plot data object

# legend_image = "04_mh_qq_plots_overview_stats/aa_excluded/legend_publication.png" #provide legend as picture because the creation with ggplot is not possible
# GET PHENOTYPES TO LOAD --------------------------------------------------

loci = fread(paste0(path_locus_definition, locusdef_file), dec = ",")
phenotypes = as.character(str_split(loci$phenotypes_in_region, " \\| ", simplify = T))
phenotypes = unique(phenotypes)
phenotypes = subset(phenotypes, phenotypes != "")
phenotypes = str_sort(phenotypes)
 


# ADDITIONAL INFORMATION FOR LABELLING ------------------------------------
loci[,cand_genes := cand_genes]
loci[,sex_ia := sex_ia]
loci[,novelty := novelty]



# CREATE PLOT DATA OBJECT -------------------------------------------------
if(collect_plot_data) {
  
  message("\n--------------------------\n")
  message("Collecting Plot Data...\n")
  
  #load meta results
  n.cores = min(length(phenotypes), max.cores)
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  tmp_out <-
    clusterCall(my.cluster, function(x)
      .libPaths(x), .libPaths())
  doParallel::registerDoParallel(cl = my.cluster)
  data = foreach (pheno = phenotypes, .packages = "data.table") %dopar% {
    pheno_file = list.files(path_data, pattern = pheno, full.names = T)
    pheno_data = fread(pheno_file, nThread = 1)
    
    #qc filtering
    pheno_data = pheno_data[pFEM < 0.01 &
                              nWeightedMAF > maf_filter &
                              nWeightedInfoScore > info_filter &
                              I2 < i2_filter  &
                              numberStudies >= n_studies_filter &
                              numberLargeStudies >= n_large_studies_filter,
                            c("markerID", "chrom", "pos", "pFEM")]
    setnames(pheno_data, "pFEM", pheno, skip_absent = T)
    pheno_data
    
  }
  parallel::stopCluster(cl = my.cluster)
  
  #bind data to plot object
  plot_data = data[[1]]
  
  if (length(data) > 1) {
    for (i in 2:length(data)) {
      plot_data = merge(plot_data,
                        data[[i]],
                        all = T,
                        by = c("markerID", "chrom", "pos"))
    }
  }
  
  #save plot object
  save(plot_data, file = paste0(path_mh_plots, data_name))
  
}


# PREPARE DATA ------------------------------------------------------------

message("\n--------------------------\n")
message("Plotting...\n")

# TODO select data
load(paste0(path_mh_plots,data_name))

#add cummulative position for x-axis labelling
plot_data <- plot_data %>% 
  group_by(chrom) %>% 
  # Compute chromosome size
  summarise(chrlength = max(pos)) %>%  
  # Calculate cumulative position of each chromosome
  mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
  select(-chrlength) %>%
  # Temporarily add the cumulative length of each chromosome to the initial 
  # dataset 
  left_join(plot_data, ., by=c("chrom"="chrom")) %>%
  # Sort by chr then position 
  arrange(chrom, pos) %>%
  # Add the position to the cumulative chromosome length to get the position of 
  # this probe relative to all other probes
  mutate(rel_pos = pos + cumulativechrlength)




# PLOT COMBINED -----------------------------------------------------------

# prepare data by converting to tall format
id_cols = c("markerID", "chrom", "pos", "cumulativechrlength", "rel_pos")
measure_cols = names(plot_data)[str_detect(names(plot_data), "_all")]

keep_cols = c(id_cols, measure_cols)
all_data = plot_data[,..keep_cols]

all_data = melt(all_data, id.vars = id_cols, measure.vars = measure_cols)
all_data = as.data.table(all_data)
all_data[,value := -log10(value)]
all_data[,variable := str_remove_all(variable, "_all")]
all_data[,variable := str_to_lower(variable)]
setnames(all_data, c("variable", "value"), c("trait", "pval"), skip_absent = T)
all_data = all_data[!is.na(pval)]

#filter p-value
maxp = ifelse(filter_max, max_p, max(all_data$pval))
if(filter_max){
  all_data[pval>max_p, pval := max_p]
}


#create labels for index variants
labels = loci[str_detect(phenotype, "_all")] #select only SNPs with all phenotype
labels[,trait := str_remove_all(phenotype, "_all")]
labels = merge(all_data, labels, by = c("markerID", "trait"), all = F)
labels[,label := paste0("locus ", region, cand_genes)]
for(r in dependent_regions){
  labels[region == r, label := str_glue("({labels[region == r, label]})")]
}


labels_sex_new = labels[sex_ia & novelty]
labels_nosex_new = labels[!sex_ia & novelty]
labels_sex_known = labels[sex_ia & !novelty]
labels_nosex_known = labels[!sex_ia & !novelty]




# positions for chrom labels
axis_df <- all_data %>% 
  group_by(chrom) %>% 
  summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)

#define the colors (alternating greys for chromosomes) and trait for gw sig hits
all_data[pval<= -log10(5 * 10 ^ -8) & chrom %% 2 == 0, color := "grey"]
all_data[pval<= -log10(5 * 10 ^ -8) & chrom %% 2 != 0, color := "grey30"]
all_data[pval> -log10(5 * 10 ^ -8), color := trait]


# plot
plot_combined <- ggplot() +
  geom_point(data = all_data[pval <= -log10(5 * 10 ^ -8)],
             aes(x = rel_pos, y = pval,
                 col = color),
             size = 0.25) +
  geom_point(data = all_data[pval > -log10(5 * 10 ^ -8)],
             aes(x = rel_pos, y = pval,
                 col = color),
             size = 1.25) +
  scale_color_manual(
    values = col_dir,
    breaks = c("cinnamon",
               "fish",
               "lemon",
               "orange",
               "pineapple",
               "score"),
    drop = F,
    name = "trait"
  ) + 
scale_x_continuous(
  labels = axis_df$chrom,
  breaks = axis_df$chr_center,
  expand = expansion(mult = 0.01)
) +
  scale_y_continuous(limits = c(2, maxp + 2 ), #add 2 for label positioning
                     expand = expansion(mult = c(0.02, 0.02))) +
  geom_hline(
    yintercept = -log10(5e-8),
    color = "red",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  labs(x = "", y = bquote(atop("OVERALL", '-log'[10] * '(p)'))) +
  theme_classic() +
  theme(legend.position = "none",
    axis.title.x = element_blank(),
    plot.margin = margin(
      b = 0,
      l = 10,
      t = 10,
      r = 5
    )) + 
  geom_label_repel(
    data = labels_nosex_new,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = c(1.9e8,rep(0,4)),
    max.overlaps = Inf,
    size = 3
  ) + 
  geom_label_repel(
    data = labels_sex_new,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    fontface = "bold",
    max.overlaps = Inf,
    size = 3
  ) +
  geom_text_repel(
    data = labels_nosex_known,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = c(0,0,1),
    nudge_x = c(-1.9e8,-2e8, 1e8),
    max.overlaps = Inf,
    size = 3
  ) + 
  geom_text_repel(
    data = labels_sex_known,
    aes(x = rel_pos, y = pval, label = label),
    fontface = "bold",
    min.segment.length = 0,
    nudge_y = 0,
    nudge_x = 1e8,
    max.overlaps = Inf,
    size = 3
  )


# PLOT FEMALES ------------------------------------------------------------

# prepare data by converting to tall format
id_cols = c("markerID", "chrom", "pos", "cumulativechrlength", "rel_pos")
measure_cols = names(plot_data)[str_detect(names(plot_data), "_female")]

keep_cols = c(id_cols, measure_cols)
female_data = plot_data[,..keep_cols]

female_data = melt(female_data, id.vars = id_cols, measure.vars = measure_cols)
female_data = as.data.table(female_data)
female_data[,value := -log10(value)]
female_data[,variable := str_remove_all(variable, "_female")]
female_data[,variable := str_to_lower(variable)]
setnames(female_data, c("variable", "value"), c("trait", "pval"), skip_absent = T)
female_data = female_data[!is.na(pval)]

#filter p-value
maxp = ifelse(filter_max, max_p, max(female_data$pval))
if(filter_max){
  female_data[pval>max_p, pval := max_p]
}


#create labels for index variants
labels = loci[str_detect(phenotype, "_female")] #select only SNPs with female top pheno
labels[,trait := str_remove_all(phenotype, "_female")]
labels = merge(female_data, labels, by = c("markerID", "trait"), all = F)
labels[,label := paste0("locus ", region, cand_genes)]
for(r in dependent_regions){
  labels[region == r, label := str_glue("({labels[region == r, label]})")]
}

labels_sex_new = labels[sex_ia & novelty]
labels_nosex_new = labels[!sex_ia & novelty]
labels_sex_known = labels[sex_ia & !novelty]
labels_nosex_known = labels[!sex_ia & !novelty]




# positions for chrom labels
axis_df <- female_data %>% 
  group_by(chrom) %>% 
  summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)

#define the colors (alternating greys for chromosomes) and trait for gw sig hits
female_data[pval<= -log10(5 * 10 ^ -8) & chrom %% 2 == 0, color := "grey"]
female_data[pval<= -log10(5 * 10 ^ -8) & chrom %% 2 != 0, color := "grey30"]
female_data[pval> -log10(5 * 10 ^ -8), color := trait]


# plot
plot_females <- ggplot() +
  geom_point(data = female_data[pval <= -log10(5 * 10 ^ -8)],
             aes(x = rel_pos, y = pval,
                 col = color),
             size = 0.25) +
  geom_point(data = female_data[pval > -log10(5 * 10 ^ -8)],
             aes(x = rel_pos, y = pval,
                 col = color),
             size = 1.25) +
  scale_color_manual(
    values = col_dir,
    breaks = c("cinnamon",
               "fish",
               "lemon",
               "orange",
               "pineapple",
               "score"),
    drop = F,
    name = "trait"
  ) + 
  scale_x_continuous(
    labels = axis_df$chrom,
    breaks = axis_df$chr_center,
    expand = expansion(mult = 0.01)
  ) +
  scale_y_continuous(limits = c(2, maxp + 2 ), #add 2 for label positioning
                     expand = expansion(mult = c(0.02, 0.02))) +
  geom_hline(
    yintercept = -log10(5e-8),
    color = "red",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  labs(x = "", y = bquote(atop("FEMALES",'-log'[10] * '(p)'))) +
  theme_classic() +
  theme(legend.position = "none",
    axis.title.x = element_blank(),
    plot.margin = margin(
      b = 0,
      l = 10,
      t = 10,
      r = 5
    )) + 
  geom_label_repel(
    data = labels_nosex_new,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    max.overlaps = Inf,
    size = 3
  ) + 
  geom_label_repel(
    data = labels_sex_new,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 4,
    nudge_x = c(-1e8,0),
    fontface = "bold",
    max.overlaps = Inf,
    size = 3
  ) +
  geom_text_repel(
    data = labels_nosex_known,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    max.overlaps = Inf,
    size = 3
  ) + 
  geom_text_repel(
    data = labels_sex_known,
    aes(x = rel_pos, y = pval, label = label),
    fontface = "bold",
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    max.overlaps = Inf,
    size = 3
  )




# PLOT MALES --------------------------------------------------------------

# prepare data by converting to tall format
id_cols = c("markerID", "chrom", "pos", "cumulativechrlength", "rel_pos")
measure_cols = names(plot_data)[str_detect(names(plot_data), "_male")]

keep_cols = c(id_cols, measure_cols)
male_data = plot_data[,..keep_cols]

male_data = melt(male_data, id.vars = id_cols, measure.vars = measure_cols)
male_data = as.data.table(male_data)
male_data[,value := -log10(value)]
male_data[,variable := str_remove_all(variable, "_male")]
male_data[,variable := str_to_lower(variable)]
setnames(male_data, c("variable", "value"), c("trait", "pval"), skip_absent = T)
male_data = male_data[!is.na(pval)]

#filter p-value
maxp = ifelse(filter_max, max_p, max(male_data$pval))
if(filter_max){
  male_data[pval>max_p, pval := max_p]
}


#create labels for index variants
labels = loci[str_detect(phenotype, "_male")] #select only SNPs with female top pheno
labels[,trait := str_remove_all(phenotype, "_male")]
labels = merge(male_data, labels, by = c("markerID", "trait"), all = F)
labels[,label := paste0("locus ", region, cand_genes)]
for(r in dependent_regions){
  labels[region == r, label := str_glue("({labels[region == r, label]})")]
}

labels_sex_new = labels[sex_ia & novelty]
labels_nosex_new = labels[!sex_ia & novelty]
labels_sex_known = labels[sex_ia & !novelty]
labels_nosex_known = labels[!sex_ia & !novelty]




# positions for chrom labels
axis_df <- male_data %>% 
  group_by(chrom) %>% 
  summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)

#define the colors (alternating greys for chromosomes) and trait for gw sig hits
male_data[pval<= -log10(5 * 10 ^ -8) & chrom %% 2 == 0, color := "grey"]
male_data[pval<= -log10(5 * 10 ^ -8) & chrom %% 2 != 0, color := "grey30"]
male_data[pval> -log10(5 * 10 ^ -8), color := trait]


# plot
plot_males <- ggplot() +
  geom_point(data = male_data[pval <= -log10(5 * 10 ^ -8)],
             aes(x = rel_pos, y = pval,
                 col = color),
             size = 0.25) +
  geom_point(data = male_data[pval > -log10(5 * 10 ^ -8)],
             aes(x = rel_pos, y = pval,
                 col = color),
             size = 1.25) +
  scale_color_manual(
    values = col_dir,
    breaks = c("cinnamon",
               "fish",
               "lemon",
               "orange",
               "pineapple",
               "score"),
    drop = F,
    name = "trait"
  ) +
  scale_x_continuous(
    labels = axis_df$chrom,
    breaks = axis_df$chr_center,
    expand = expansion(mult = 0.01),
    position = "top"
  ) +
  scale_y_reverse(limits = c(maxp + 2, 2 ), #add 2 for label positioning
                     expand = expansion(mult = c(0.02, 0.02))) +
  geom_hline(
    yintercept = -log10(5e-8),
    color = "red",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  labs(x = "", y = bquote(atop("MALES",'-log'[10] * '(p)'))) +
  theme_classic() +
  theme(legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = ggplot2::element_blank(),
    plot.margin = margin(
      b = 10,
      l = 10,
      t = 0,
      r = 5
    )) + 
  geom_label_repel(
    data = labels_nosex_new,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    max.overlaps = Inf,
    size = 3
  ) + 
  geom_label_repel(
    data = labels_sex_new,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    fontface = "bold",
    max.overlaps = Inf,
    size = 3
  ) +
  geom_text_repel(
    data = labels_nosex_known,
    aes(x = rel_pos, y = pval, label = label),
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    max.overlaps = Inf,
    size = 3
  ) + 
  geom_text_repel(
    data = labels_sex_known,
    aes(x = rel_pos, y = pval, label = label),
    fontface = "bold",
    min.segment.length = 0,
    nudge_y = 1,
    nudge_x = 0,
    max.overlaps = Inf,
    size = 3
  )




# CREATE DUMMY PLOT FOR LEGEND THAT CONTAINS ALL VALUES -------------------

dummy_dat = rbind(all_data, female_data, male_data)
dummy_legend_plot <- ggplot() +
  geom_point(data = dummy_dat[pval <= -log10(5 * 10 ^ -8)],
             aes(x = rel_pos, y = pval,
                 col = trait),
             size = 3) +
  scale_color_manual(
    values = col_dir,
    limits = c("cinnamon",
               "fish",
               "lemon",
               "orange",
               "pineapple",
               "score"),
    drop = F,
    name = "traits"
  ) + theme(legend.position="bottom") +
  guides(colour = guide_legend(nrow = 1)) #force one row


legend <- g_legend(dummy_legend_plot)

# ARRANGE FINAL PLOT ------------------------------------------------------

miami_plot = grid.arrange(plot_females, plot_males, nrow = 2)
final_plot = grid.arrange(plot_combined, miami_plot, legend, nrow = 3, heights = c(1,1,0.2))
ggsave(paste0(path_mh_plots, image_name), final_plot, height = 16, width = 22, units = "cm")

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
