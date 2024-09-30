#' ---
#' title: "Locus Definition"
#' subtitle: "Define loci for all future analysis"
#' author: "Katrin Horn, Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' Definition of the genome wide significant loci
#'
#' # Initialize ####
#' ***
#' Update server & set working directory
# pipeline_name: 05a_locus_definition.R

rm(list = ls())
time0 = Sys.time()

source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server


setwd(projectpath)

#' # Parameters for QC ####

p_filter = 5 * 10 ^ -8
maf_filter = 0.01
info_filter = 0.8 
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1 # only kept to be compatible with analysis pipeline, no filtering results from this in the final analysis

#' # Analysis ####
#' ***
#' ## get genome-wide sig. SNPs ####
#'
#' 1. Step: get genome-wide significant SNPs of all phenotypes and settings

message("\n--------------------------\n")
message("Loading data...\n")

myFiles = list.files(path = path_data, pattern = ".gz")

result.1 = foreach(f = myFiles) %do% {
  input_fn = paste0(path_data, f)
  dat = fread(input_fn, nThread = 30)
  
  #add a phenotype column
  phenotype = str_match(f, "GWASMA_([^_]+_[^_]+)_\\d")[2]
  dat[, phenotype := rep(phenotype, nrow(dat))]
  
  # make a column that counts the large studies (by subtracting ARIC from the total count)
  # aric_present = (!is.na(dat$p.ARIC_AFR)) + (!is.na(dat$p.ARIC_EUR))
  # dat[,numberLargeStudies := dat$numberStudies - aric_present]
  
  #filter for genome-wide significance, MAF and number of studies and info (no longer needed when working on filtered data)
  dat = dat[pFEM < p_filter &
              nWeightedMAF > maf_filter &
              nWeightedInfoScore > info_filter &
              I2 < i2_filter  &
              numberStudies >= n_studies_filter &
              numberLargeStudies >= n_large_studies_filter, ]
  
  #' add a position column (the script assumes a continuous position here over all
  #' chromosomes). Therefore it is assumed that each chromosome
  #' is 250.000.000 bp long (length of chromosome 1). This should make sure that loci
  #' on different chromosomes are not mistaken as neighbors.
  dat[, position := chrom * 250000000 + pos]
  dat
}
result.1 = rbindlist(result.1)
#table(result.1[, phenotype])

# helper function to get minimal distance between significant SNPs --> do we have to collapse or not?
getSmallestDist = function(x) {
  #if there is only one significant SNP we want to beak the while loop
  if (length(x) <= 1) {
    return(500000)
  }
  
  y = c(x[2:length(x)], max(x) + 1000000)
  return(min(y - x))
}

#' ## Shrink SNP list by 500 kb range ###
#'
#' 2. Step: shrink list to best SNPs in phenotype and setting (loop: take best SNP and remove SNPs in +/- 500 kb range until no SNPs left to remove from range of a better SNP)
subsets = unique(result.1[, phenotype])

result.2 = foreach(s = subsets) %do% {
  #s = subsets[1]
  subdata = copy(result.1)
  subdata = subdata[phenotype == s,]
  setkey(subdata, position)
  subdata[, keep := NA]
  subdata[, NR_SNPs := as.numeric(NA)]
  
  smallestDist = getSmallestDist(subdata[, position])
  while (smallestDist < 500000) {
    minP = min(subdata[is.na(keep), pFEM])
    subdata[minP == pFEM, keep := T]
    myPos = subdata[minP == pFEM, position]
    stopifnot(length(myPos) == 1)
    #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
    myFilt = (subdata[, position] < (myPos - 500000)) |
      (subdata[, position] > (myPos + 500000)) |
      subdata[, keep]
    myFilt[is.na(myFilt)] = FALSE
    subdata = subdata[myFilt == TRUE,]
    #NR_SNPs sums the filtered SNPs for a locus
    subdata[minP == pFEM, NR_SNPs := sum(myFilt == F)]
    smallestDist = getSmallestDist(subdata[, position])
  }
  
  #stopifnot(sum(is.na(subdata[,keep])) <= 1)
  subdata[is.na(keep), NR_SNPs := 0]
  subdata[is.na(keep), keep := TRUE]
  subdata
}
result.2 = rbindlist(result.2)

#add regions
result.2[, region_start := position - 500000]
result.2[, region_end := position + 500000]
#result.2

#' ## Shrink SNP list overall ###
#'
#' 3. Step: again shrink list to best SNPs of either phenotype
#'
#'  Same procedure as step 3. See how regions collapse over phenotypes.
subdata = copy(result.2)
setkey(subdata, position)

#assign new region number to all lines of data.table
subdata[, region := as.numeric(NA)]
subdata[1, region := 1]
done = foreach(l = c(2:nrow(subdata))) %do% {
  if (subdata[l, region_start] <= subdata[l - 1, region_end]) {
    subdata[l, region := subdata[l - 1, region]]
  } else {
    subdata[l, region := subdata[l - 1, region] + 1]
  }
}

#keep best SNPs per region
allRegions = unique(subdata[, region])

result.4 = foreach(r = allRegions) %do% {
  subdata.2 = copy(subdata)
  subdata.2 = subdata.2[region == r,]
  minP = min(subdata.2[, "pFEM"])
  subdata.2[pFEM == minP, region_start := min(subdata.2[, region_start])]
  subdata.2[pFEM == minP, region_end := max(subdata.2[, region_end])]
  subdata.2[pFEM == minP, NR_SNPs := sum(subdata.2[, NR_SNPs], na.rm = T) + nrow(subdata.2) - 1]
  subdata.2[pFEM == minP, phenotypes_in_region := paste(unique(subdata.2[, phenotype]), collapse = " | ")]
  subdata.2[pFEM == minP,]
}
result.4 = rbindlist(result.4)

#correct the region coordinates back to chromosome wise coordinates
result.4[, "region_start"] = result.4[, region_start] - result.4[, chrom] * 250000000
result.4[, "region_end"] = result.4[, region_end] - result.4[, chrom] * 250000000


cols2keep = c(
  "region",
  "region_start",
  "region_end",
  "markerID",
  "phenotype", 
  "chrom",
  "pos",
  "aa",
  "ea",
  "nWeightedMAF",
  "nWeightedInfoScore",
  "I2",
  "betaFEM",
  "direction",
  "seFEM",
  "pFEM",
  "NR_SNPs",
  "numberStudies",
  "numberLargeStudies",
  "totalN",
  "phenotypes_in_region"
)

# add the single study characteristics
study_columns = merge(c("p.", "beta.", "se.", "n.", "maf.", "infoscore."),
                      c("LIFE_EUR", "Rhineland_EUR","ARIC_EUR", "CHRIS_EUR")) #TODO change studies if necessary
study_columns = paste0(study_columns$x, study_columns$y)

cols2keep = c(cols2keep, study_columns)

colsOut = setdiff(colnames(result.4), cols2keep)
result.4[, get("colsOut") := NULL]
setcolorder(result.4, cols2keep)
#result.4

message("\n--------------------------\n")
message(" Writing locus definiton file...\n")

write.table(
  result.4,
  file = paste0(path_locus_definition, "locus_definition.csv"), 
  col.names = T,
  row.names = F,
  quote = F,
  sep = ";",
  dec = ","
)


message("\n--------------------------\n")
message("Finished locus definition.\n")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "mins"), 3), " minutes")
