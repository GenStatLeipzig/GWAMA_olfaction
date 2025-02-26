# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-08-07
#
# Script Description: Perform a power analysis of the different odour traits
#
#
# Notes:
#
# pipeline_name: T9_power_calculation.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R")
library(genpwr)
library(scales)
setwd(projectpath)


# VARIABLES ---------------------------------------------------------------


# maf = 0.01
maf_list = c(0.01, 0.05, 0.1, 0.2)
betas = seq(0.01, 0.5, 0.01)
or = exp(betas)
a = 5E-8
n = 18895


f_case_rate = "additional_information/case_rates.xlsx"

# DATA PREPARATION --------------------------------------------------------


# calculate case rate based on all studies
d = read_excel(f_case_rate) %>% as.data.table()

d = d[subgroup == "all"]
d[,correct_total := correct_LIFE + correct_CHRIS + correct_ARIC + correct_RHINELAND]
d[,incorrect_total := incorrect_LIFE + incorrect_CHRIS + incorrect_ARIC + incorrect_RHINELAND]
d[,case_rate := incorrect_total/(correct_total + incorrect_total)]
# setorder(d, -case_rate)

case_rates = d[["case_rate"]]
names(case_rates) = d[["odour"]]



# POWER CALCULATION -------------------------------------------------------
for (maf in maf_list) {
  f_out = paste0("03_power_plots/power_plot_maf_", maf, ".pdf")
  res = foreach(odds_rat = or) %do% {
    res = foreach(p = names(case_rates)) %do% {
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
      pow[, pheno := p]
      pow[, OR := odds_rat]
    }
    res = rbindlist(res)
  }
  res = rbindlist(res)


  # VISUALIZATION -----------------------------------------------------------
  setorder(res, -Case.Rate)
  res[, pheno := str_remove_all(pheno, "_all")]
  res[, label := str_glue("{res$pheno} ({round(res$Case.Rate,2)})")]
  col = viridis(12, direction = -1)
  setnames(res, paste0("Power_at_Alpha_", a), "Power")


  ggplot(data = res, aes(x = OR, y = Power, colour = label)) +
    geom_line(linewidth = 1.5) +
    scale_color_manual(breaks = unique(res$label), values = col) +
    labs(title = str_glue("Alpha: {a}, MAF: {maf}, N: {n}")) +
    guides(color = guide_legend(title = "Trait (Case rate)")) +
    scale_x_continuous(
      trans = log_trans(),
      breaks = exp(c(0, 0.1, 0.2, 0.3, 0.4, 0.5)),
      labels = label_math(e^.x, format = log)
    ) +
    ylab("Power") +
    ylim(c(0, 1))

  ggsave(f_out)
}

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

# sessionInfo()
message("\nTOTAL TIME : ", round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
