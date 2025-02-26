# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-02
#
# Script Description: Plot correlation between the different odours based on LIFE-Adult data
#
#
# Notes:
#
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------
p_out = "17_smell_correlation/"
directory = p_out
if (!dir.exists(directory)) {
  dir.create(directory)
}


life_data = fread("phenoFile.txt") # contains individual level results from odour detection. File will not be provided.

# translate odours from German to English
translation = c(
  "Orange" = "orange",
  "Schuhleder" = "leather",
  "Zimt" = "cinnamon",
  "Pfefferminz" = "peppermint",
  "Banane" = "banana",
  "Zitrone" = "lemon",
  "Lakritz" = "liquorice",
  "Kaffee" = "coffee",
  "Gewuerznelke" = "cloves",
  "Ananas" = "pineapple",
  "Rose" = "rose",
  "Fisch" = "fish",
  "SCORE" = "SCORE"
)
phenotypes = names(life_data) %>%
  str_split_i(., "_", 1) %>%
  unique() %>%
  .[!str_detect(., "SCORE|FID|IID")]

# exclude score at it is continuous
col_all_de = names(life_data)[str_detect(names(life_data), pattern = "_all") &
  !str_detect(names(life_data), pattern = "SCORE")]
col_male_de = names(life_data)[str_detect(names(life_data), pattern = "_males") &
  !str_detect(names(life_data), pattern = "SCORE")]
col_female_de = names(life_data)[str_detect(names(life_data), pattern = "_females") &
  !str_detect(names(life_data), pattern = "SCORE")]

col_all = purrr::pmap(list(string = col_all_de, pattern = phenotypes, replacement = translation[phenotypes]), str_replace_all) %>% unlist()
col_male = purrr::pmap(list(string = col_male_de, pattern = phenotypes, replacement = translation[phenotypes]), str_replace_all) %>% unlist()
col_female = purrr::pmap(list(string = col_female_de, pattern = phenotypes, replacement = translation[phenotypes]), str_replace_all) %>% unlist()

setnames(
  life_data,
  c(col_all_de, col_male_de, col_female_de),
  c(col_all, col_male, col_female)
)

col_all = sort(col_all)
col_male = sort(col_male)
col_female = sort(col_female)

# CALCULATE CORRELATION MATRICES ------------------------------------------

#  phi would be the correct metric to use. However Pearson's r^2 is a reasonable good approximation considering our case numbers
cor_all = cor(life_data[, ..col_all])
cor_female = cor(life_data[, ..col_female], use = "pairwise.complete.obs")
cor_male = cor(life_data[, ..col_male], use = "pairwise.complete.obs")

dimnames(cor_all) = dimnames(cor_all) %>% purrr::map(., ~ str_replace_all(.x, "_", " "))
dimnames(cor_female) = dimnames(cor_female) %>% purrr::map(., ~ str_replace_all(.x, "_", " "))
dimnames(cor_male) = dimnames(cor_male) %>% purrr::map(., ~ str_replace_all(.x, "_", " "))

# VISUALIZATION -----------------------------------------------------------

pdf(paste0(p_out, "smell_correlation_all_publication.pdf"))
pheatmap(cor_all, display_numbers = TRUE, main = "Correlation of smell (r) - overall", cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

pdf(paste0(p_out, "smell_correlation_female_publication.pdf"))
pheatmap(cor_female, display_numbers = TRUE, main = "Correlation of smell (r) - female", cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

pdf(paste0(p_out, "smell_correlation_male_publication.pdf"))
pheatmap(cor_male, display_numbers = TRUE, main = "Correlation of smell (r) - male", cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
