# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-04-29
#
# Script Description: Creates LDSR command for single study statistics
#
#
# Notes:
# pipeline_name: R11b_create_ldsr_commands.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R") # TODO check server
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

files = list.files("19_ldsr_single_studies/sum_stats/", pattern = ".sumstats.gz")
traits = str_match(files, "^(.*)\\.pipeline")[, 2]

template = "helper_scripts/ldsr_template.txt"
script = "19_ldsr_single_studies/ldsr_command.sh"

directory = "19_ldsr_single_studies/heritability/"
if (!dir.exists(directory)) {
  dir.create(directory)
}

# CREATE LDSC COMMAND -----------------------------------------------------

write("", script)
for (trait in traits) {
  infile = files[str_detect(files, trait)]
  infile = paste0("19_ldsr_single_studies/sum_stats/", infile)
  outfile = paste0("19_ldsr_single_studies/heritability/", trait, "_h2")

  ldsr_command = readLines(template)
  ldsr_command = str_replace_all(ldsr_command, "FILELIST", infile)
  ldsr_command = str_replace_all(ldsr_command, "OUT", outfile)

  write(ldsr_command, script, append = TRUE)
  write("\n", script, append = TRUE)
}
system(paste0("chmod +x ", script))

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
