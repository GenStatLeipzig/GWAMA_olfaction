# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-04-29
#
# Script Description: Creates LDSR command
#
#
# Notes:
# pipeline_name: R02b_create_ldsr_commands.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

files = list.files("14_heritability/sum_stats_formatted/", pattern = ".sumstats.gz")
traits = str_match(files, "GWASMA_(.*)_\\d+-")[, 2]

template = "helper_scripts/ldsr_template.txt"
script = "14_heritability/ldsr_command.sh"

directory = "14_heritability/heritability/"
if (!dir.exists(directory)) {
  dir.create(directory)
}

# Create command------------------------------------------------------------

write("", script)
for (trait in traits) {
  infile = files[str_detect(files, trait)]
  infile = paste0("14_heritability/sum_stats_formatted/", infile)
  outfile = paste0("14_heritability/heritability/", trait, "_h2")

  ldsr_command = readLines(template)
  ldsr_command = str_replace_all(ldsr_command, "FILELIST", infile)
  ldsr_command = str_replace_all(ldsr_command, "OUT", outfile)

  write(ldsr_command, script, append = T)
  write("\n", script, append = T)
}
system(paste0("chmod +x ", script))

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
