# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-04-26
#
# Script Description: Prepares summary statistics of individual studies
# for LD score regression
# RsID is attatched based on hg38 reference by matching via chromosome and bp.
# SNPs with mismatching alleles are filtered out later, when using the preparation
# step for LDSR and merging with the LD reference with consideration of the alleles.
#
#
# Notes:
# munge script has to be executed manually.
#
# pipeline_name: R11a_sumstats_preparation.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

path_data = "" # overwrite default data path with single study summary stats

chroms_to_load = c(1:22)
ld_sum_stats_folder = "19_ldsr_single_studies/sum_stats_unformatted/"
ld_sum_stats_folder_formatted = "19_ldsr_single_studies/sum_stats/"
file_list = list.files(path_data)

cols_to_keep = c(
  "ID",
  "codedAll",
  "noncodedAll",
  "n",
  "p",
  "maf",
  "infoscore",
  "beta",
  "se",
  "chr",
  "pos"
)
new_names = c(
  "snpid",
  "A1",
  "A2",
  "N",
  "P-value",
  "maf",
  "info",
  "beta",
  "se",
  "chrom",
  "pos"
)

munge_template = "helper_scripts/munge_command_template.txt"
munge_command_file = "19_ldsr_single_studies/munge_command.sh"

maf_filter = 0.01
info_filter = 0.8 # note that the default value of LDSR is 0.9

# LOAD CHROMOSOME REFERENCE DATA ------------------------------------------
# for hg38

message("\n--------------------------\n")
message("Loading data...\n")

n.cores = min(40, length(chroms_to_load))
my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
snp_annotation = foreach(
  chrom = chroms_to_load,
  .packages = c("data.table", "stringr"),
  .export = c("path_snp_annotation")
) %dopar% {
  if (chrom == 23) chrom = "X"

  annotation = fread(str_glue("{path_snp_annotation}homo_sapiens-chr{chrom}.vcf.gz"),
    nThread = 1
  )
  annotation = annotation[, c("ID", "#CHROM", "POS")]
  annotation
}
parallel::stopCluster(cl = my.cluster)

snp_annotation = rbindlist(snp_annotation)
setnames(snp_annotation, "#CHROM", "CHROM", skip_absent = T)
setkey(snp_annotation, CHROM, POS)




# PROCESS SINGLE STUDIES SUM STATS FILES ----------------------------------

# add rsID, select columns and change their names to match with the tool
# do not use parallel as copying the annotation data takes to much time

result = foreach(f = file_list, .packages = c("data.table")) %do% {
  data = fread(paste0(path_data, f), nThread = 30)
  data[, maf := ifelse(eaf > 0.5, 1 - eaf, eaf)]

  # qc filter the data to save some time
  data = data[maf > maf_filter &
    infoscore > info_filter]

  # attatch rsID via chrom and bp position
  data[, chrPosId := paste(chr, pos, sep = ":")]
  setkey(data, chr, pos)
  data = merge(data, snp_annotation[, c("ID", "CHROM", "POS")], by.x = c("chr", "pos"), by.y = c("CHROM", "POS"), all = FALSE)

  # select and change some columns
  data = data[, ..cols_to_keep]
  setnames(data, new_names)

  directory = ld_sum_stats_folder
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
  directory = ld_sum_stats_folder_formatted
  if (!dir.exists(directory)) {
    dir.create(directory)
  }

  # output
  fwrite(
    data,
    paste0(ld_sum_stats_folder, f),
    sep = " ",
    col.names = TRUE,
    row.names = FALSE,
    nThread = 30
  )
}

# CREATE FILE FOR MUNGE STEP ----------------------------------------------

command = readLines(munge_template)
write("", munge_command_file, append = F)

for (f in file_list) {
  mungecom = copy(command)
  mungecom = str_replace_all(mungecom, "STATS", paste0(ld_sum_stats_folder, f))
  mungecom = str_replace_all(mungecom, "OUT", paste0(ld_sum_stats_folder_formatted, str_remove(f, ".gz")))

  write(mungecom, munge_command_file, append = T)
  write("\n", munge_command_file, append = T)
}

system("chmod +x 19_ldsr_single_studies/munge_command.sh")

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")
message("TODO: Manually activate conda envionment and run the munge_command.sh script")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
