# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-04-26
#
# Script Description: Prepares summary statistics for coffee and tea intake for bivariate LDSC
#
#
# Notes:
# no info score provided, SNPs are filtered for MAF > 0.01, no qc has to be performed
# No infoscore should not cause problems, because we filter to HapMap3 in the munge step
#
# munge script has to be executed manually
#
# GWAS summary statistics for coffee odour identification need to be manually copied to directory
#
# pipeline_name: R12a_ldsc_preparation.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

path_data = "18_ldsr_coffee/sum_stats_ukbb/" # files downloaded from https://yanglab.westlake.edu.cn/pub_data.html
ld_sum_stats_folder = "18_ldsr_coffee/sum_stats_raw/"
ld_sum_stats_folder_formatted = "18_ldsr_coffee/sum_stats_formatted/"
file_list = list.files(path_data, pattern = "\\.gz$")

# data is in hg19. As ldsc matches via rsID that should be identical, no lifting is required
cols_to_keep = c(
  "SNP",
  "A1",
  "A2",
  "N",
  "P",
  "maf",
  "b",
  "se"
)
new_names = c(
  "snpid",
  "A1",
  "A2",
  "N",
  "P-value",
  "maf",
  "beta",
  "se"
)

munge_template = "helper_scripts/munge_command_template.txt"
munge_command_file = "18_ldsr_coffee/munge_command.sh"

# PROCESS GWAMA SUM STATS FILES -------------------------------------------

purrr::walk(c(ld_sum_stats_folder, ld_sum_stats_folder_formatted), ~ if (!dir.exists(.x)) {
	dir.create(.x, recursive = TRUE)
})

result = foreach(f = file_list, .packages = c("data.table")) %do% {
  data = fread(paste0(path_data, f), nThread = 30)
	data[, maf := ifelse(freq>0.5, 1-freq, freq)]

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

system("chmod +x 18_ldsr_coffee/munge_command.sh")

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")
message("TODO: Manually activate conda envionment and run the munge_command.sh script")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
