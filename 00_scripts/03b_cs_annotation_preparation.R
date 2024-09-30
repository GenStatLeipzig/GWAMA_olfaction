# preparation of the annotation input file
# data is lifted from hg38 to hg19 coordinates because pipeline is currently 
# based on that
# SNPs are taken until a sumProb of >99% is reached
# creates an annotation file per CS instead of all CS combined to handle stack overflow
# errors in the pipeline
# pipeline_name: 10b_cs_annotation_preparation_single_regions.R

####INIT####
rm(list = ls())

source("00_scripts/00_SourceFile_smelling_meta.R")
source("helper_scripts/liftOverHg38TOHg19.R")
setwd(projectpath)

path_credible_set_unlifted = "08_credible_set_analysis/"
path_output = "08_credible_set_analysis/annotation_input/"


filelist = list.files(path_credible_set_unlifted, pattern = "region_\\d")


####PROCESSING####
process_file = function(file) {
  credible_set_data = fread(paste0(path_credible_set_unlifted, file))
  
  pos = min(which(credible_set_data[,SumProb]>0.99))
  credible_set_data = credible_set_data[c(1:pos), ]
  
  lifted_data = data.table("chr" = credible_set_data$Chr, "pos" = credible_set_data$bp, "snps" = credible_set_data$SNP)
  
  lifted_data = liftOVerHg38TOHg19(
    lifted_data,
    try_unlifted_via_ensembl = F,
    return_as_data_table = T,
    path_chainfile = path_chainfile
  )
  
  annotation_data = data.table("snp" = lifted_data$snps, "chr_hg19" = lifted_data$chr_hg19, "pos_hg19" = lifted_data$pos_hg19)
  annotation_data = merge(annotation_data, credible_set_data[,c("SNP", "PostProb", "SumProb", "refA", "freq", "n", "p","info")], by.x = "snp", by.y = "SNP", all.x=T, all.y=F, )
  annotation_data[, minor_allele_freq := sapply(annotation_data$freq, calculate_maf_from_freq)]
  annotation_data[,region := str_match(file, "region_(\\d+)_")[2]]
  annotation_data[,phenotype := str_match(file, "region_\\d+_([^_]+_[^_]+)_\\d+")[2]]
  
  #colnames(credible_set_data)[c(1, 4, 5)] = c("ID", "CHROM", "GENPOS")
  credible_set_name = str_replace(file, ".txt", "")
  annotation_data[, credible_set := credible_set_name]
  
  
  for (row in 1:nrow(annotation_data)){
    alleles = str_match_all(annotation_data[row, snp], "^\\d+:\\d+:(.+):(.+)$")
    allele1 = alleles[[1]][,2]
    allele2 = alleles[[1]][,3]
    
    if (allele1==annotation_data[row, refA]){
      annotation_data[row, a2:=allele2]
    } else if (allele2 == annotation_data[row, refA]){
      annotation_data[row, a2:=allele1]
    } else {
      message("Mismatching alleles")
    }
  }
  
  setnames(annotation_data, "refA", "a1", skip_absent = T)
  output = annotation_data[,c("snp", "chr_hg19", "pos_hg19", "a1", "a2", "minor_allele_freq", "credible_set","region","phenotype", "n", "p", "info", "PostProb", "SumProb")]
  setnames(output, "p", "P-value", skip_absent = T)
  
  fout = str_glue("snp_annotation_{credible_set_name}.txt")
  write.table(
    output,
    paste0(path_output, fout),
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
  )
  
  return(output)
  
}

annotation_data = foreach(
  file = filelist,
  .packages = c(
    "data.table",
    "GenomicRanges",
    "rtracklayer",
    "reshape2",
    "stringr"
  )
) %do% {
  process_file(file)
}

annotation_data = rbindlist(annotation_data)
write.table(
  annotation_data,
  paste0(path_output, "snp_annotation_all_credible_sets.txt"),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

print("Finished Preparing Annotation File")