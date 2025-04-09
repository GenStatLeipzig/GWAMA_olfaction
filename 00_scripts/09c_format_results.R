# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-12-09
#
# Script Description: Convert .log file into csv format
#
#
# Notes:
# pipeline_name: R06c_format_results.R
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
library(openxlsx)

# VARIABLES ---------------------------------------------------------------

start_row = "p1"
end_row = "Analysis finished"

p_analysis = "15_ldsr_diseases/"
p_ukbb = "UKBB_phenotypes/"
p_out = paste0(p_analysis, "genetic_correlation/")


f_phenoslct = paste0("additional_information/", "ldsc_selected_phenotypes.xlsx")
f_phenoslct_f = paste0("additional_information/", "ldsc_selected_phenotypes_females.xlsx")
f_phenoslct_m = paste0("additional_information/", "ldsc_selected_phenotypes_males.xlsx")
f_out = paste0(p_out, "genetic_correlation_multiple_phenotypes.xlsx")

new_col_names = c(
	"p1" = "phenotype1",
	"p2" = "phenotype2",
	"sex" = "pheno2.sex",
	"rg" = "genetic correlation",
	"se" = "se genetic correlation",
	"p" = "p",
	"p_fdr" = "p_fdr",
	"h2_obs" = "pheno2.h2",
	"h2_obs_se" = "pheno2.h2_se",
	"h2_int" = "pheno2.h2_int",
	"h2_int_se" = "pheno2.h2_int_se",
	"phenotype" = "pan-UKBB identifier",
	"description" = "description",
	"source" = "source",
	"ldsc_sumstat_dropbox" = "ldsc_sumstat_dropbox"
)

# PROCESSING --------------------------------------------------------------

fl = list.files(p_out, pattern = "\\.log$", full.names = TRUE)

d = foreach(f_in = fl) %do% {
  d = readLines(f_in)
  start = which(str_detect(d, start_row))
  end = which(str_detect(d, end_row)) - 2 # -2 due to the empty line above
  d = d[start:end]

  writeLines(d, "tmp_table.csv")
  d = fread("tmp_table.csv")
  file.remove("tmp_table.csv")
    d
}
d = rbindlist(d)



# FORMATTING --------------------------------------------------------------

d[,p1 := str_match(p1, "GWASMA_(.*)_\\d+")[,2]]
d[, p2 := str_split_i(p2, "__", 1)]
d[, p2 := str_split_i(p2, "//", 2)]
d[1, p2 := "score_all"]



# COMBINE WITH STUDY INFORMATION ------------------------------------------

d[,sex := str_split_i(p1, "_", 2)]
d[sex == "all", sex := "both_sexes"]

d_ukbb = read_excel(f_phenoslct) %>% as.data.table()
d_ukbb_f = read_excel(f_phenoslct_f) %>% as.data.table()
d_ukbb_m = read_excel(f_phenoslct_m) %>% as.data.table()

# format descriptors to match on LDSC results
format_descriptors = function(d_tmp){
	d_tmp <- d_tmp[!is.na(description)]
	
matching_description = ifelse(
		str_starts(d_tmp$description, "Non-cancer illness code"),
		str_split_i(d_tmp$description, ": ", 2),
		d_tmp$description
) %>%
  tolower() %>%
  str_replace_all(., " ", "_")
	d_tmp[, matching_description := matching_description]
	stopifnot(all(d_tmp$matching_description %in% d$p2))

	d_tmp = d_tmp[, .(phenotype, description, source, sex, ldsc_sumstat_dropbox, matching_description)]
	return(d_tmp)
}

d_ukbb = foreach(d_tmp = list(d_ukbb, d_ukbb_f, d_ukbb_m)) %do% {
	format_descriptors(d_tmp)
}
d_ukbb = rbindlist(d_ukbb)


d = merge(d, d_ukbb, by.x = c("p2", "sex"), by.y = c("matching_description", "sex"), all.x = TRUE, sort = FALSE)
cols = names(d)
cols[1] = "p1"
cols[2] = "p2"
cols[3] = "sex"
d = d[, ..cols]
setorder(d, -rg, na.last = TRUE)

d[,p_fdr := p.adjust(d$p, method = "fdr")]
d = d[,.SD, .SDcols = names(new_col_names)]
setnames(d, names(new_col_names), new_col_names)

WriteXLS::WriteXLS(
  d,
  str_replace(f_out, "\\.xlsx$", "_supplement.xlsx"),
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
