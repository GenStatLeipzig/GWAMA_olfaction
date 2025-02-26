library(data.table)

df <- as.data.frame(fread("coloc_data_cropped_regions_qc_filtered.csv.gz"))

# for each region, extract start and end position
res <- c()
for (region in unique(df$region)){
    sub <- df[which(df$region==region),]
    print(length(unique(sub$chrom)))
    #print(length(unique(sub$phenotype)))
    res <- rbind(res, c(region, unique(sub$chrom), min(sub$pos), max(sub$pos)))
}

# might be multiple phenotypes
res <- as.data.frame(res)
write.table(res, file="coloc_region.txt", row.names=F, col.names=c("region", "chr", "start", "end"), quote=F, sep="\t")
