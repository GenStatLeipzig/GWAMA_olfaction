# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2023-11-24
#
# Script Description: Performing a meta regression with MR-MEGA
#
#
# Notes:
# Alleles with more than 10 nucleotides are filtered out
# This method only provides p-values and no betas (given betas correspond to primary
# components and not the effect of the alleles)
#
# pipeline name: R10a_mrmega_regression.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------
mrmega_path = ".../MRMEGA/" # path to MR-MEGA executable
path_reformatted_data_pipeline = "" # path to results from single studies in a common format
path_analysis = "13_heterogeneity/mrmega/" # analysis folder
input_files = list.files(path_reformatted_data_pipeline)

# limit analysis to phenotypes with genome-wide significant signals
loci = fread(
  paste0(
    path_locus_definition,
    "locus_definition_rsID.csv"
  ),
  dec = ","
)
phenos = loci[, phenotypes_in_region] %>%
  str_split(., " \\| ") %>%
  unlist() %>%
  unique()
input_files = input_files[str_detect(input_files, str_flatten(phenos, collapse = "|"))]

convert_files = TRUE # set to true when the data has to be converted into the mr_mega input format
n.cores = min(length(input_files), 32)
pc = 1 # number of primary components (must be < n_studies-2, i.e. for 5 cohorts it can be at most 2)

qc_filter = FALSE # could filter out variants that have sufficient n-weighted MAF or infoscore
maf_filter = 0.01
info_filter = 0.8

# CONVERT FORMAT ----------------------------------------------------------

if (convert_files == T) {
  message("\n--------------------------\n")
  message("Converting Input Files...\n")


  my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
  tmp_out = clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
  doParallel::registerDoParallel(cl = my.cluster)
  result = foreach(f = input_files, .packages = c("data.table", "stringr")) %dopar% {
    data = fread(paste0(path_reformatted_data_pipeline, f), nThread = 1)

    if (qc_filter) {
      data[, maf := ifelse(eaf > 0.5, 1 - eaf, eaf)]
      data = data[maf > maf_filter & infoscore > info_filter]
    }

    if (str_detect(f, "SCORE")) {
      # formatting for quantitatvie traits with effect ans SE instead of OR
      setnames(data,
        c("markerID", "codedAll", "noncodedAll", "beta", "se", "eaf", "n", "chr", "pos"),
        c("MARKERNAME", "EA", "NEA", "BETA", "SE", "EAF", "N", "CHROMOSOME", "POSITION"),
        skip_absent = F
      )
    } else {
      # formatting of qualitative traits
      setnames(data,
        c("markerID", "codedAll", "noncodedAll", "eaf", "n", "chr", "pos"),
        c("MARKERNAME", "EA", "NEA", "EAF", "N", "CHROMOSOME", "POSITION"),
        skip_absent = F
      )
      data[, OR := exp(beta)]
      data[, OR_95L := exp(beta - se * 1.96)]
      data[, OR_95U := exp(beta + se * 1.96)]
    }

    # some alleles observed in ARIC had a confusing structure where
    # the alternative allele was long (caused errors when reading files into MR-MEGA)
    # filter these SNPs out
    data = data[str_length(EA) <= 10 & str_length(NEA) <= 10, ]

    # save
    f_out = str_replace(f, "pipeline", "mrmega")

    directory = paste0(path_analysis, "/input")
    if (!dir.exists(directory)) {
      dir.create(directory)
    }

    fwrite(
      data,
      paste0(path_analysis, "input/", f_out),
      row.names = FALSE,
      quote = F,
      sep = "\t",
      nThread = 1
    )
  }
  parallel::stopCluster(cl = my.cluster)
}


# META REGRESSION ---------------------------------------------------------
message("\n--------------------------\n")
message("Perfomring MR-MEGA Regression...\n")

directory = paste0(path_analysis, "output/")
if (!dir.exists(directory)) {
  dir.create(directory)
}

phenotypes = phenos

setwd(mrmega_path) # change to program directory due to issues with the call (./MR-MEGA can not be integrated into path)
n.cores = min(length(phenotypes), 32)
my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out = clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
result = foreach(phenotype = phenotypes, .packages = c("data.table", "stringr"), .export = "pc") %dopar% {
  # create list of input files
  phenotype_files = list.files(paste0(projectpath, path_analysis, "input/"), pattern = phenotype, full.names = T)
  mrmega_infile = paste0(projectpath, path_analysis, "MR_MEGA_", phenotype, ".in")
  write(phenotype_files, file = mrmega_infile)

  mrmega_outfile = paste0(projectpath, path_analysis, "output/", "mrmega_", phenotype)
  # run mrmega
  if (str_detect(phenotype, "SCORE")) {
    # add qt flag for quantitative traits
    mrmega_command = str_glue("./MR-MEGA --pc {pc} -i {mrmega_infile} -o {mrmega_outfile} --qt")
  } else {
    mrmega_command = str_glue("./MR-MEGA --pc {pc} -i {mrmega_infile} -o {mrmega_outfile}")
  }
  system(mrmega_command)
  file.remove(mrmega_infile)
}
parallel::stopCluster(cl = my.cluster)

# PLOTTING AND P-VALUE CORRECTION -----------------------------------------

message("\n--------------------------\n")
message("Plotting MH and Q-Q-plots...\n")
setwd(projectpath)

directory = paste0(path_analysis, "plots/")
path_mrmega_plots = directory
if (!dir.exists(directory)) {
  dir.create(directory)
}

file_list = list.files(paste0(path_analysis, "output"), full.names = T, pattern = "\\.result$")

n.cores = min(32, length(file_list))
my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out = clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
result = foreach(f = file_list, .packages = "stringr", .export = "mrmega_path") %dopar% {
  # in the script exists a lower bound for the p-values of  p-values>1e-14 due to library
  # reasons. This will be corrected by the following step to p-values>1e-325
  correction_call = str_glue("R --slave --vanilla --args input={f} out={f} < {mrmega_path}fixP.r")
  system(correction_call)

  # plotting
  phenotype = str_match(f, "mrmega_(.+_.+)\\.result$")[2]
  mh_out = paste0(path_mrmega_plots, phenotype, "_mh_plot.png")
  qq_out = paste0(path_mrmega_plots, phenotype, "_qq_plot.png")

  mh_call = str_glue("R --slave --vanilla --args input={f} out={mh_out} < {mrmega_path}manh.r")
  system(mh_call)

  qq_call = str_glue("R --slave --vanilla --args input={f} out={qq_out} < {mrmega_path}qq.r")
  system(qq_call)
}
parallel::stopCluster(cl = my.cluster)

# FILE COMPRESSION --------------------------------------------------------

message("\n--------------------------\n")
message("Compressing Files...\n")

file_list = list.files(paste0(path_analysis, "output"), full.names = T, pattern = "\\.result$")

n.cores = min(32, length(file_list))
my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out = clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
result = foreach(f = file_list, .packages = "stringr") %dopar% {
  gzip_call = str_glue("gzip {f}")
  system(gzip_call)
}
parallel::stopCluster(cl = my.cluster)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")


# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
