# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2025-02-07
#
# Script Description: Add gene names to mehtylation probes.
#
#
# Notes:
#
# pipeline_name: R09_add_methylation_genes.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
setwd(projectpath)

# VARIABLES ---------------------------------------------------------------

setwd("COLOC") # enter coloc path here
f_qtl_coloc = "2025_02_06_olfaction_qtl_colocs.xlsx"
f_methylation_annot = "humanmethylation450_15017482_v1-2.csv" # Illumina methylation 450k reference from manifest. Needs to be downloaded from Illumina website.
f_out = "2025_02_06_olfaction_qtl_colocs_annotated_methylation.xlsx"

ref = 1 # which reference to use (1 or 2)

# annotation --------------------------------------------------------------

d_qtl = read_excel(f_qtl_coloc) %>% as.data.table(.)


d_annot = fread(f_methylation_annot, skip = "IlmnID", fill = TRUE)
d_annot = d_annot[str_starts(IlmnID, "cg")]


d = merge(d_qtl, d_annot[, .(IlmnID, UCSC_RefGene_Name)], by.x = "probe", by.y = "IlmnID", all.x = TRUE, sort = FALSE)
d[probe.gene == "NA", probe.gene := UCSC_RefGene_Name]
d[, UCSC_RefGene_Name := NULL]
setorder(d, "region_index", "qtl", "phenotype", "probe.gene")

WriteXLS(
  d,
  str_replace(f_out, ".xlsx", "_supplement.xlsx"),
  AdjWidth = TRUE,
  AutoFilter = TRUE,
  BoldHeaderRow = TRUE,
  AllText = TRUE
)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
