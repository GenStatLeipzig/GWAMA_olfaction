# first iterate over all results files and extract significant results

th = 0.8
res = NULL
nr.snps = 50

files = list.files(path = "results", full.names = TRUE, recursive = TRUE)

for (file in files) {
  cat("collecing file", file, "of", length(files), "\n")
  df = read.delim(file, header = T)
  i = which(df$PP.H4.abf >= th & df$nsnps >= nr.snps)
  if (length(i) > 0) {
    res = rbind(res, df[i, ])
  }
}

# add probe info
map = read.delim("epi_mapping.txt", header = T)
m = match(res$probe, map$probe_id)
res = cbind(res, probe.chr = map$chr[m], probe.position = map$probe_bp[m], probe.gene = map$gene[m])


write.table(res, file = "all_significant_colocs_coloc_abf.txt", row.names = F, col.names = T, sep = "\t", quote = F)
