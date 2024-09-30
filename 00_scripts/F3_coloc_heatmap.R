# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
# 
# Date: 2024-01-11
#
# Script Description: Visualization of the coloc results in form of a heatmap
# Is formatted for publication with annotated description of loci
#
#
# Notes: 
# TODO before executing the scrip manually add the results for significance and
# higher sex effect from the sex.ia test output to the coloc results
# pipeline_name: 11c_coloc_visualization_publication.R
#

# INIT --------------------------------------------------------------------

rm(list=ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
suppressPackageStartupMessages(library(pheatmap))
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------

# coloc_file = "11_sex_interaction/coloc_sex_ia.csv"
coloc_file = "11_sex_interaction/coloc_sex_ia.csv" #TODO check file

outfile = "11_sex_interaction/coloc_publication.pdf"


# PUBLICATION READY PLOT --------------------------------------------------

coloc = fread(coloc_file, dec = ",")

# remove region 7 as it is not independent
coloc = coloc[-c(7)]

rnames = transpose(str_split(coloc$locus, " "))[[2]]
plot_data = expand.grid(locus = rnames, coloc_pp = c("H0", "H1", "H2", "H3", "H4"))
plot_data$pp = c(coloc$PP.H0.abf, coloc$PP.H1.abf, coloc$PP.H2.abf, coloc$PP.H3.abf, coloc$PP.H4.abf)
number_cols = ifelse(plot_data$pp<0.5, "white", "black")

# change colour of region labels according to results of sex IA test
x_col = coloc$sex.higherEffect
x_col = replace(x_col, which(x_col == "none"), "black")
x_col = replace(x_col, which(x_col == "female"), "#B2182B")
x_col = replace(x_col, which(x_col == "male"), '#2166AC')
x_col[6] = "black" #TODO the code above colours all sex IA, also the ones that are not FDR sig, remove color manually for them

# change font type to make sexIA annotation better visible
x_face = coloc$sex.higherEffect
x_face = replace(x_face, which(x_face == "none"), "italic")
x_face = replace(x_face, which(x_face == "female"), "bold")
x_face = replace(x_face, which(x_face == "male"), 'bold')
x_face[6] = "italic"

ggplot(plot_data, aes(locus, coloc_pp, fill = pp)) +
  geom_tile(show.legend = F, color = "black") +
  scale_fill_viridis(discrete = F) +
  geom_text(aes(label = round(pp, 2)), col = number_cols) + 
  xlab("locus") + 
  ylab("COLOC hypothesis") +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme(axis.ticks = element_blank(),
        axis.title.x = element_text(face="bold", size = 12),
        axis.title.y = element_text(face="bold", size = 12),
        axis.text.x = element_text(color = x_col, size = 11, face=x_face),
        axis.text.y = element_text(color = "black", size = 11))

ggsave(outfile, height = 7, width = 17, units = "cm")

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
