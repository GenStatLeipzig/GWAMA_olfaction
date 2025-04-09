# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-08-07
#
# Script Description: Perform a power analysis of the different odour traits.
# Power is calculated for different population sizes.
#
#
# Notes:
#
#

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
library(genpwr)
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------
f_case_rate = "additional_information/case_rates.xlsx"

maf_list = c(0.01, 0.05, 0.1, 0.2)
betas = seq(0.01, 0.5, 0.01)
or = exp(betas)
a = 5E-8

case_numbers = c(
  Meta_all = 18895,
  Meta_male = 8757,
  Meta_female = 10138,
  Gisladottir_2020 = 11326,
  Dong_2017 = 1979 + 6582,
  McRae_2013 = 187
)


# DATA PREPARATION --------------------------------------------------------

d = read_excel(f_case_rate) %>% as.data.table()

d = d[subgroup == "all"]
d[, correct_total := correct_LIFE + correct_CHRIS + correct_ARIC + correct_RHINELAND]
d[, incorrect_total := incorrect_LIFE + incorrect_CHRIS + incorrect_ARIC + incorrect_RHINELAND]
d[, case_rate := incorrect_total / (correct_total + incorrect_total)]
# setorder(d, -case_rate)

case_rates = d[["case_rate"]]
names(case_rates) = d[["odour"]]

directory = "03_power_plots/"
if (!dir.exists(directory)) {
  dir.create(directory)
}

# POWER CALCULATION AND VISUALIZATION -------------------------------------

for (maf in maf_list) {
  n.cores = min(30, length(names(case_rates)))
  my.cluster = parallel::makeCluster(n.cores, type = "FORK")
  # tmp_out <- clusterCall(my.cluster, function(x) .libPaths(x), .libPaths())
  doParallel::registerDoParallel(cl = my.cluster)
  foreach(p = names(case_rates)) %dopar% {
    ## data collection ####
    phenotype = str_remove(p, "_all$")
    f_out = paste0(directory, "power_plot_samplesize_maf_", maf, "_", phenotype, ".pdf")

    res = foreach(odds_rat = or) %do% {
      res = foreach(n_i = seq_along(case_numbers)) %do% {
        n = case_numbers[n_i]
        pow = genpwr.calc(
          calc = "power",
          model = "logistic",
          N = n,
          MAF = maf,
          Alpha = a,
          Case.Rate = case_rates[p],
          OR = odds_rat,
          Test.Model = "Additive",
          True.Model = "Additive"
        )

        pow = as.data.table(pow)
        pow[, study := paste0(str_replace_all(names(case_numbers)[n_i], "_", " "), ": ", n)]
        pow[, n := n]
        pow[, OR := odds_rat]
      }
      res = rbindlist(res)
    }
    res = rbindlist(res)


    # visualization ####
    setorder(res, -n)
    setnames(res, paste0("Power_at_Alpha_", a), "Power")
    col = nord::nord("aurora", n = length(case_numbers))

    ggplot(data = res, aes(x = OR, y = Power, colour = study)) +
      geom_line(linewidth = 1.5) +
      scale_color_manual(breaks = unique(res$study), values = col) +
      labs(title = str_glue("Alpha: {a}, MAF: {maf}, Case rate: {round(case_rates[p], 2)} ({phenotype})")) +
      guides(color = guide_legend(title = "Sample size")) +
      scale_x_continuous(
        trans = log_trans(),
        breaks = exp(c(0, 0.1, 0.2, 0.3, 0.4, 0.5)),
        labels = label_math(e^.x, format = log)
      ) +
      ylab("Power") +
      ylim(c(0, 1))

    ggsave(f_out, width = 15, height = 10, units = "cm")
  }
  parallel::stopCluster(cl = my.cluster)
}

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
