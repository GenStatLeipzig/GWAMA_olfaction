# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-11-28
#
# Script Description:
# Plot histograms of I^2 distribution for each trait.
# Collect data in tables
#
#
#
# Notes:
# pipeline_name: R01a_heterogeneity.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)
library(pheatmap)

# VARIABLES ---------------------------------------------------------------

max.cores = 40

maf_filter = 0.01
info_filter = 0.8
n_studies_filter = 2

path_out = "13_heterogeneity/"
directory = path_out
if (!dir.exists(directory)) {
  dir.create(directory)
}


# FUNCTIONS ---------------------------------------------------------------

get_heterogeneity_snps = function(pheno, plot = FALSE) {
  d = snp_data[[pheno]]
  d[, qc_pass := nWeightedMAF > maf_filter &
    nWeightedInfoScore > info_filter &
    numberStudies >= n_studies_filter]
  if (plot) {
    ggplot() +
      geom_histogram(mapping = aes(I2, fill = qc_pass), data = d) +
      ggtitle(str_glue("I2 distribution"), pheno)
    ggsave(str_glue("{path_out}histogramm_I2_{pheno}.png"))
  }


  snp_count = list(
    "I2 QC >= 85" = d[I2 >= 85 & qc_pass == TRUE, .N] / d[qc_pass == TRUE, .N]
  )
  return(snp_count)
}

# LOAD DATA ---------------------------------------------------------------

fl = list.files(path_data, "\\.gz$")
pheno_names = str_match(fl, "GWASMA_(.*)_\\d")[, 2]

n.cores = min(length(fl), max.cores)
my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out = clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
snp_data = foreach(infile = fl, .packages = loadedNamespaces()) %dopar% {
  fread(paste0(path_data, infile), nThread = 1)
}
parallel::stopCluster(cl = my.cluster)
names(snp_data) = pheno_names


# COLLECT HETEROGENEITY MEASURES ------------------------------------------

i2_numbers = lapply(pheno_names, get_heterogeneity_snps)
rownames = names(i2_numbers[[1]])
names(i2_numbers) = pheno_names
i2_numbers = lapply(i2_numbers, unlist)
setDT(i2_numbers)

i2_numbers[, row_names := rownames]


fwrite(
  i2_numbers,
  paste0(path_out, "fraction_of_heterogeneous_snps.csv"),
  row.names = FALSE,
  quote = FALSE,
  sep = ";",
  dec = ","
)


# make plot for publication
plot_data = i2_numbers[row_names == "I2 QC >= 85"]
select_cols = names(plot_data) %>% str_detect(., "_all")
plot_data = plot_data[, ..select_cols]
names(plot_data) = names(plot_data) %>%
  str_to_lower() %>%
  str_remove(., "_all")
plot_data = plot_data * 100 # convert to percent

# sorting columns
sorting = as.vector(plot_data) %>% unlist()
sorting = sort(sorting)
sorting = names(sorting)
plot_data = plot_data[, ..sorting]

display_numbers = round(plot_data, 4) %>%
  paste0(., "%") %>%
  matrix(., nrow = 1)
pdf(paste0(path_out, "heterogeneity_snp_fraction_publication.pdf"), width = 8, height = 2.3)
pheatmap(
  plot_data,
  color = viridis(10, begin = 0, end = 0.25),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = display_numbers,
  show_rownames = FALSE,
  number_color = "white"
  # number_format = "%.4%"
)
dev.off()

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
