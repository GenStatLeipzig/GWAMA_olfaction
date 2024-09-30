# HEADER ------------------------------------------------------------------
#
# Author: Franz Förster
#
# Date: 2024-07-16
#
# Script Description: eQTL-Coloc
#
#
# Notes: Extension of the original analysis. Instead of only investigating genes with
# reported eQTLs in LD to index SNP in their corresponding tissues, aditionally all genes in 
# a 250kb window around the index variant are also analysed. For all genes, all tissues
# will be analysed
#
# Performs COLOC between observed odour associations and eQTLs to check whether
# observed signals can further influence the expression of other genes.
# For each locus genes are selected when the index SNP is in LD (r^2 > 0.3) with
# eQTLs from GTEx or the gene is in a 250kb window around the index variant.
# All available tissues will be tested.
# TODO some Ensemble IDs have to be checked manually (marked in gene table)

# General workflow:
# 1. Create list of genes regulated by variants in LD with index variants + near genes
# 2. Select Data of olfactory associations and gene expression associations
# 3. Apply QC criteria, filter for respective gene, find overlap between oflactory
# and GE snp set
# 4. perform COLOC

# pipeline_name: 18_comprehensive_eqtl_coloc.R

# INIT --------------------------------------------------------------------
rm(list = ls())
time0 = Sys.time()
source("00_scripts/00_SourceFile_smelling_meta.R") #TODO check server
setwd(projectpath)
library(coloc)

# VARIABLES ---------------------------------------------------------------

f_locusdef = "locus_definition.csv" #TODO check locus definition file
p_annotation = "08_credible_set_analysis/results/topliste_tabdelim/" #TODO check annotation folder of eqtl data

f_gene_id_conversion = "additional_information/hgnc_071624.txt" #custom hgnc download to convert gene symbol to ensembl ID (downloaded on 16.07.24, 14:47)
custom_assignment = c("ZBED9" = "ENSG00000232040", "ZNRD1" = "ENSG00000066379") #manual assignment for genes that could not be matched uniquely to EnsemblID

f_gtx_anno = "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"

max.cores = 40

loci = fread(paste0(path_locus_definition, f_locusdef), dec = ",")
gene_id_conv = fread(f_gene_id_conversion)

log_file = "10_eqtl_coloc/log_progress.txt"

maf_filter = 0.01
info_filter = 0.8
i2_filter = 85
n_studies_filter = 2
n_large_studies_filter = 1

tissues_complete = list.files(p_gtx, "GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_.*\\.allpairs.txt.gz")
tissues_complete = str_match(tissues_complete, "all_associations_(.*)\\.allpairs.txt.gz$")[,2]


# FUNCTIONS ---------------------------------------------------------------

list_eqtl_genes = function(r) {
  region_id = loci[r, region]
  index_snp = loci[r, markerID]
  
  
  #select proximate genes
  f_anno = list.files(p_annotation,
                      pattern = str_glue("proximate_genes_.*_region_{region_id}_.*\\.txt$"),
                      full.names = T)
  anno = fread(f_anno)
  genes1 = anno[markername == index_snp]$genename
  genes1 = unique(genes1)
  
  # extract genes linked though eQTLs
  f_anno = list.files(
    p_annotation,
    pattern = str_glue("eqtlinfo_.*_region_{region_id}_.*\\.txt$"),
    full.names = T
  )
  anno = fread(f_anno)
  
  genes = anno[snps == index_snp &
                 `cistrans` == "cis" &
                 study == "GTEx_V8"] #TODO select wanted studies
  genes = genes$genesymbol
  genes = unique(genes)
  genes = data.table(genesymbol = c(genes1, genes))
  genes = unique(genes)
  genes = genes[genesymbol != ""]
  
  if(length(genes) == 0){
    return()
  }
  
  # Ensembl ID for each gene (start by direct match of gene symbol than search for
  # aliases if none found)
  for (g in 1:nrow(genes)) {
    ensblID = gene_id_conv[`Approved symbol` == genes[g, genesymbol], `Ensembl gene ID`]
    if (genes[g, genesymbol] %in% names(custom_assignment)) {
      ensblID = custom_assignment[genes[g, genesymbol]]
    }
    if (length(ensblID) == 0) {
      ensblID = gene_id_conv[str_detect(`Previous symbols` , genes[g, genesymbol]), 
                             `Ensembl gene ID`]
    }
    if (length(ensblID) == 0) {
      ensblID = gene_id_conv[str_detect(`Alias symbols` , genes[g, genesymbol]), 
                             `Ensembl gene ID`]
    }
    if (length(ensblID) == 0) {
      ensblID = "NA"
    }
    if (length(ensblID) > 1) {
      ensblID = paste(ensblID, collapse = "|")
      genes[g, ensemblID := ensblID]
      genes[g, manual_decision := T] #mark if a manual descion about the ID is needed
    } else {
      genes[g, ensemblID := ensblID]
      genes[g, manual_decision := F]
    }
  }
  
  # add locus data
  genes[, region := region_id]
  genes[, snp := index_snp]
  genes[, chrom := loci[r, chrom]]
  genes[, start := loci[r, region_start]]
  genes[, end := loci[r, region_end]]
  genes[, beta_olf := loci[r, betaFEM]]
  genes[, se_olf := loci[r, seFEM]]
  genes[, p_olf := loci[r, pFEM]]
  genes[, pheno := loci[r, phenotype]]
  
  setcolorder(
    genes,
    c(
      "region",
      "snp",
      "chrom",
      "start",
      "end",
      "beta_olf",
      "se_olf",
      "p_olf",
      "pheno",
      "genesymbol",
      "ensemblID",
      "manual_decision"
    )
  )
  return(genes)
}

perform_eqtl_coloc = function(r){
  
  d_smell = smelling_data[[gene_table[r, pheno]]]
  d_smell = d_smell[(chrom == gene_table[r, chrom]) &
                      (pos >= gene_table[r, start]) &
                      (pos <= gene_table[r, end]) &
                      (numberStudies >= n_studies_filter) &
                      (numberLargeStudies >= n_large_studies_filter) &
                      (I2 < i2_filter) &
                      (nWeightedInfoScore > info_filter) &
                      (nWeightedMAF > maf_filter)]
  
  
  d_gtx = tissue_data[[gene_table[r, tissue]]]
  d_gtx = d_gtx[maf > maf_filter &
                  chrom == gene_table[r, chrom] &
                  ensemblID == gene_table[r, ensemblID]]
  
  #remove duplicates
  dupl = d_gtx$pos[duplicated(d_gtx$pos)]
  d_gtx = d_gtx[!(pos %in% dupl)]
  
  # collect common snps by matching of position and alleles
  common_snps = merge(d_smell[,c("markerID", "ea", "aa", "pos")], d_gtx[,c("variant_id", "a1", "a2", "pos")], by = "pos", all = F)
  common_snps[,allele_match := F]
  common_snps[a1 == ea & a2 == aa, allele_match := T]
  common_snps[a1 == aa & a2 == ea, allele_match := T]
  common_snps = common_snps[allele_match == T]
  
  if(nrow(common_snps) == 0){
    result = data.table(
      e_snp = NA,
      beta_exp = NA,
      se_exp = NA,
      p_exp = NA,
      n_exp = NA,
      a1_e = NA,
      a2_2 = NA,
      trait1 = NA,
      trait2 = NA,
      nsnps = 0,
      PP.H0.abf = NA,
      PP.H1.abf = NA,
      PP.H2.abf = NA,
      PP.H3.abf = NA,
      PP.H4.abf = NA
    )
    
    sink(log_file, append = T)
    cat("Done row ",r, "\n")
    sink(file = NULL) #close the log file
    
    return(result)
  }
  
  d_smell = d_smell[markerID %in% common_snps$markerID]
  d_gtx = d_gtx[variant_id %in% common_snps$variant_id]
  
  setorder(d_smell, pos)
  setorder(d_gtx, pos)
  
  d_smell[, chrPosID := paste(chrom, pos, sep = ":")]
  d_gtx[, chrPosID := paste(chrom, pos, sep = ":")]
  
  # prepare data to coloc format
  gwas1 = list(
    pvalues = d_smell$pFEM,
    N = max(d_smell$totalN),
    MAF = d_smell$nWeightedMAF,
    beta = d_smell$betaFEM,
    varbeta = (d_smell$seFEM) ^ 2,
    type = ifelse(str_detect(gene_table[r,pheno], "(?i)score"), "quant", "cc"),
    snp = d_smell$chrPosID
  )
  
  gwas2 = list(
    pvalues = d_gtx$pval_nominal,
    N = max(d_gtx$n_samples),
    MAF = d_gtx$maf,
    beta = d_gtx$slope,
    varbeta = (d_gtx$slope_se) ^ 2,
    type = "quant",
    snp = d_gtx$chrPosID
  )
  
  
  # check_dataset(gwas1)
  # check_dataset(gwas2)
  coloc_result = coloc.abf(gwas1, gwas2)
  
  #collect data about eQTL SNP
  e_snp = d_gtx[pval_nominal == min(d_gtx$pval_nominal), variant_id]
  if(length(e_snp)>1){
    e_snp = e_snp[1]
  }
  e_snp_data = d_gtx[variant_id == e_snp]
  
  result1 = data.table(
    e_snp = e_snp,
    beta_exp = e_snp_data$slope,
    se_exp = e_snp_data$slope_se,
    p_exp = min(d_gtx$pval_nominal),
    n_exp = e_snp_data$n_samples,
    a1_e = e_snp_data$a1,
    a2_2 = e_snp_data$a2,
    trait1 = gene_table[r, pheno],
    trait2 = paste("GE", gene_table[r, genesymbol])
  )
  result2 = as.data.table(as.list(coloc_result$summary))
  result = cbind(result1, result2)
  
  sink(log_file, append = T)
  cat("Done row ",r, "\n")
  sink(file = NULL) #close the log file
  
  return(result)
}

# COLLECT GENE LIST -------------------------------------------------------

# for each locus select genes where eQTLs are in LD with the respective index variant
# or are in physical proximity to index variant
# only eQTLs reported by GTEx are considered (LIFE data requires publication)

message("\n--------------------------\n")
message("Create genelist...\n")

gene_table = foreach (r = 1:nrow(loci)) %do% {
  list_eqtl_genes(r)
}

gene_table = rbindlist(gene_table)
message(
  str_glue(
    "{length(unique(gene_table[ensemblID == 'NA']$genesymbol)) + length(unique(gene_table[ensemblID == '']$genesymbol))} genes will be discarded to to missing ID. See 'discarded_genes.txt' for details."
  )
)

writeLines(unique(gene_table[ensemblID == "NA" | ensemblID == ""]$genesymbol),
           "10_eqtl_coloc/discarded_genes.txt")

gene_table = gene_table[ensemblID != "NA"]
gene_table = gene_table[ensemblID != ""] 

fwrite(
  gene_table,
  "10_eqtl_coloc/gene_table.csv",
  row.names = FALSE,
  quote = F,
  sep = ";",
  dec = ","
)




# PREPARE DATA AND PERFORM COLOC ------------------------------------------
message("\n--------------------------\n")
message("Preparing phenotype data ...\n")
# add tissues (all tissues should be checked for all genes)
gene_table = fread("10_eqtl_coloc/gene_table.csv", dec = ",")

tg = expand.grid(unique(gene_table$genesymbol), tissues_complete)
tg = as.data.table(tg)
names(tg) = c("genesymbol", "tissue")
setkey(tg, "genesymbol")

gene_table = merge(gene_table, tg, by = "genesymbol", all = F, allow.cartesian = T)
gene_table[,tissue := as.character(tissue)]

# load smelling data
phenos = unique(gene_table$pheno)

n.cores = min(length(phenos), max.cores)
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
tmp_out <- clusterCall(my.cluster, function(x)
  .libPaths(x), .libPaths())
doParallel::registerDoParallel(cl = my.cluster)
smelling_data = foreach(pheno = phenos,
                        .packages = c("data.table", "stringr")) %dopar% {
                          f_pheno = list.files(path_data, pheno)
                          d = fread(paste0(path_data, f_pheno), nThread = 1)
                        }
parallel::stopCluster(cl = my.cluster)

names(smelling_data) = phenos

# as memory is to limited to load all gtex tissues at once we will do the following process
# for a subset of tissues that is loaded in parallel
gene_table_ori = copy(gene_table)
tissue_set_1 = tissues_complete[1:10]
tissue_set_2 = tissues_complete[11:20]
tissue_set_3 = tissues_complete[21:30]
tissue_set_4 = tissues_complete[31:40]
tissue_set_5 = tissues_complete[41:length(tissues_complete)]
i = 1

for(tissue_filter in list(tissue_set_1, tissue_set_2, tissue_set_3, tissue_set_4, tissue_set_5)){

  message("\n--------------------------\n")
  message("Prepare eQTL data of tissueset ", i, "...\n")
  
  
  gene_table = gene_table_ori[tissue %in% tissue_filter]
  

  # load GTEx data for selected tissues
  tissues = unique(gene_table$tissue)
  
  n.cores = min(length(tissues), max.cores)
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  tmp_out <- clusterCall(my.cluster, function(x)
    .libPaths(x), .libPaths())
  doParallel::registerDoParallel(cl = my.cluster)
  tissue_data = foreach (t = tissues, .packages = c("data.table", "stringr")) %dopar% {
    f_tissue = list.files(p_gtx, pattern = str_glue("{t}\\.allpairs\\.txt\\.gz$"))
    d = fread(paste0(p_gtx, f_tissue), nThread = 1)
    d[, ensemblID := str_remove(gene_id, "\\..+$")]
    
    pos_data = transpose(str_split(d$variant_id, "_"))
    d[, chrom := str_remove(pos_data[[1]], "chr")]
    d[, chrom := as.numeric(str_replace(chrom, "X", "23"))]
    d[, pos := as.numeric(pos_data[[2]])]
    d[, a1 := pos_data[[3]]]
    d[, a2 := pos_data[[4]]]
    d[, n_samples:=round((ma_count/maf)/2,1)]
  }
  parallel::stopCluster(cl = my.cluster)
  names(tissue_data) = tissues
  
  
  # PERFRORM COLOCALIZATION -------------------------------------------------
  
  # apply QC filter (MAF>1% for GTEx, the ususal QCs for smelling data)
  # select overlap of SNPs between smelling and GTEx data for the respective gene
  # perform COLOC
  # lifting is no longer needed as both datasets use hg38
  
  message("\n--------------------------\n")
  message("performing COLOC...\n")
  
  sink(log_file, append = F)
  cat("TODO: ", nrow(gene_table), " entries\n")
  sink(file = NULL) #close the log file
  
  
  # do not use parallel processing as copying the eQTL data takes to much memory
  coloc_result = foreach (r = 1:nrow(gene_table), .packages = c("data.table", "stringr", "coloc")) %do% {
    perform_eqtl_coloc(r)
  }
  
  coloc_result = rbindlist(coloc_result)
  final_output = cbind(gene_table, coloc_result)
  
  fwrite(
    final_output,
    str_glue("10_eqtl_coloc/eqtl_coloc_tissueset_{i}.csv"),
    row.names = FALSE,
    quote = F,
    sep = ";",
    dec = ","
  )
  
  i = i+1
  
  
}

# TODO read all the result files, combine them and sort them (make an interactive excel file)


# COMBINE THE RESULTS -----------------------------------------------------

tissue_files = list.files("10_eqtl_coloc/", "eqtl_coloc_tissueset", full.names = T)

result = foreach(f = tissue_files) %do% {
  d = fread(f, dec = ",")
}
result = rbindlist(result)

setorder(result, region, genesymbol, tissue)

WriteXLS(result, "10_eqtl_coloc/eqtl_coloc_combined.xlsx",AdjWidth = T, AutoFilter = T, BoldHeaderRow = T)

# END ---------------------------------------------------------------------
message("\n--------------------------\n")
message("Finished.\n")

#sessionInfo()
message("\nTOTAL TIME : " , round(difftime(Sys.time(), time0, units = "hours"), 2), " hours")
