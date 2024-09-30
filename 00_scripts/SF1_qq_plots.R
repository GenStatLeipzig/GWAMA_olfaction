# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2023-06-22
#
# Script Description: Creates QQ-Plots for supplementary
#
# Notes:
# Prints the Q-Q-plots ready for publication (multiple plots per page, only QC SNPs are plotted)
# works in sequential manner for easier plotting
# pipeline_name: 04e_qq_plots_publication.R

rm(list = ls())
source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
n.cores = 40 #cores for parallelization 

setwd(projectpath)
library(gridExtra)
library(grid)

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

ci = 0.95

reducePointsToPlot = T # reduce number of points by random selection

#### function definition

summary_analysis <- function(file) {
  phenotype_data = fread(paste0(path_data, file), nThread = 1) 
  
  # get analysis data
  phenotype = str_match(file, "GWASMA_(\\D+)_\\d")[2]
  phenotype_name = str_split(phenotype, "_", simplify = T)[1]
  analysis_name = str_split(phenotype, "_", simplify = T)[2]
  
  plot_data_qq = phenotype_data[!is.na(pFEM) &
                                      nWeightedMAF > maf_filter &
                                      nWeightedInfoScore > info_filter &
                                      I2 < i2_filter  &
                                      numberStudies >= n_studies_filter &
                                      numberLargeStudies >= n_large_studies_filter,
                                    c("markerID", "chrom", "pos", "pFEM")]
  plot_data_qq[,color := "black"]
  
  lambda = phenotype_data[!is.na(pFEM) &
                            nWeightedMAF > maf_filter &
                            nWeightedInfoScore > info_filter &
                            I2 < i2_filter  &
                            numberStudies >= n_studies_filter &
                            numberLargeStudies >= n_large_studies_filter,
                          median((betaFEM / seFEM) ^ 2, na.rm = TRUE) / qchisq(0.5, 1)]
  lambda = round(lambda, 4)
  
  
  plot_data_qq = as.data.frame(plot_data_qq)
  
  # build QQ-plot
  df = data.table::data.table(pvalues = plot_data_qq$pFEM,
                              qc = plot_data_qq$qc_snp,
                              color = plot_data_qq$color)
  data.table::setorder(df, pvalues)
  n  <- nrow(df)
  df[, observed := -log10(pvalues)]
  df[, expected := -log10(ppoints(n))]
  df[, clower   := -log10(qbeta(
    p = (1 - ci) / 2,
    shape1 = 1:n,
    shape2 = n:1
  ))]
  df[, cupper   := -log10(qbeta(
    p = (1 + ci) / 2,
    shape1 = 1:n,
    shape2 = n:1
  ))]
  
  # reduce the number of points
  if (reducePointsToPlot == TRUE) {
    random = runif(n, min=0, max=1)
    df[, random := random]
    df[observed < 6, toPlot := random>0.75]
    df[observed < 5, toPlot := random>0.90]
    df[observed < 3, toPlot := random>0.99]
    df[observed < 2, toPlot := random>0.9999]
    df[observed >= 6, toPlot := TRUE]
    df = df[toPlot == TRUE, ]
    
  } 
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  qqplot = ggplot2::ggplot(df, ggplot2::aes(x = expected, y = observed, color = qc)) +
    ggplot2::geom_ribbon(
      mapping = ggplot2::aes(
        x = expected,
        ymin = clower,
        ymax = cupper,
        color = "black"
      ),
      alpha = 0.1
    ) + ggplot2::geom_point(ggplot2::aes(expected, observed, col = "black"), size = 1) +
    ggplot2::geom_abline(intercept = 0,
                         slope = 1,
                         alpha = 0.5) +
    ggplot2::xlab(log10Pe) +
    ggplot2::ylab(log10Po) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_line(size = 0.5),
      panel.grid = ggplot2::element_blank(),
      # legend.position = c(0.2, 0.7)
      legend.position = "none"
      # panel.grid = element_line(size = 0.5, color = "grey80")
    ) +
    ggplot2::scale_color_manual(values = c("black")) +
    ggplot2::scale_shape_manual(values = c(1, 2)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(8)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    ggplot2::ggtitle(paste0(str_to_lower(phenotype_name), " ", analysis_name),
                     subtitle = bquote(lambda == .(lambda))) #unicode \u03BB
  
  
  return(qqplot)
   

}

####parallel processing####
#preparation
file_list = list.files(paste0(path_data), "\\.gz$")


my.cluster <- parallel::makeCluster(n.cores,
                                    type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
out_lib <- clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())

#parallel processing
plots = foreach (file = file_list,
                        .packages = c("data.table", "stringr", "CMplot")) %dopar% {
                          summary_analysis(file)
                        }
parallel::stopCluster(cl = my.cluster)

# save to pdf
# grid.arrange(rectGrob(), rectGrob())
multiplot = marrangeGrob(plots, ncol = 3, nrow = 5, layout_matrix = matrix(1:15, 5, 3, TRUE), top = NULL)
ggsave(paste0(path_mh_plots, "qq_plots_publication.pdf"), multiplot,  width = 210, height = 297, units = "mm")



print("Finished Q-Q-Plots")