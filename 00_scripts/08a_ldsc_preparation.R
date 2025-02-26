# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-04-26
#
# Script Description: Prepares the meta-summary statistics for LD score regression
# RsID is attached based on hg38 reference by matching via chromosome and bp.
# SNPs with mismatching alleles are filtered out later, when using the preparation
# step for LDSR and merging with the LD reference with consideration of the alleles.
#
#
# Notes:
# munge script has to be executed manually:
# ./16_LDSC/munge_command.sh
#
# pipeline_name: R02a_ldsc_preparation.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

chroms_to_load = c(1:22)
ld_sum_stats_folder = "14_heritability/sum_stats_raw/"
ld_sum_stats_folder_formatted = "14_heritability/sum_stats_formatted/"
file_list = list.files(path_data)

cols_to_keep = c(
  "ID",
  "ea",
  "aa",
  "totalN",
  "pFEM",
  "nWeightedMAF",
  "nWeightedInfoScore",
  "betaFEM",
  "seFEM",
  "chrom",
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
munge_command_file = "14_heritability/munge_command.sh"


qc_filter = TRUE
maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1


# LOAD CHROMOSOME REFERENCE DATA ------------------------------------------
# for hg38
purrr::walk(c(ld_sum_stats_folder, ld_sum_stats_folder_formatted), ~ if (!dir.exists(.x)) {
  dir.create(.x, recursive = TRUE)
})


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
  annotation
}
parallel::stopCluster(cl = my.cluster)

snp_annotation = rbindlist(snp_annotation)
# snp_annotation[, chrPosId := paste(`#CHROM`, POS, sep = ":")]
setnames(snp_annotation, "#CHROM", "CHROM", skip_absent = T)
setkey(snp_annotation, CHROM, POS)

# remove unneded columns to save disk space
snp_annotation = snp_annotation[, c("ID", "CHROM", "POS")]

# PROCESS GWAMA SUM STATS FILES -------------------------------------------
# add rsID, select columns and change their names to match with the tool

# do not use parallel as copying the annotation data takes to much time
result = foreach(f = file_list, .packages = c("data.table")) %do% {
  data = fread(paste0(path_data, f), nThread = 30)

  # qc filter the data to save some time
  if (qc_filter) {
    data = data[nWeightedMAF > maf_filter &
      nWeightedInfoScore > info_filter &
      I2 < i2_filter &
      numberStudies >= n_studies_filter &
      numberLargeStudies >= n_large_studies_filter, ]
  }

  # attach rsID via chrom and bp position
  data[, chrPosId := paste(chrom, pos, sep = ":")]
  setkey(data, chrom, pos)
  data = merge(data, snp_annotation[, c("ID", "CHROM", "POS")], by.x = c("chrom", "pos"), by.y = c("CHROM", "POS"), all = F)

  # select and change some columns
  data = data[, ..cols_to_keep]
  setnames(data, new_names)

  # output
  fwrite(
    data,
    paste0(ld_sum_stats_folder, f),
    nThread = 1,
    sep = " ",
    col.names = T,
    row.names = F
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

system("chmod +x 14_heritability/munge_command.sh")

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")
message("TODO: Manually activate conda envionment and run the munge_command.sh script")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
