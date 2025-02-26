args <- commandArgs(trailingOnly=T)

esi <- args[1]
chr <- as.numeric(args[2])
start <- as.numeric(args[3])
end <- as.numeric(args[4])
out <- args[5]

df <- read.delim(esi, header=F)
i <- which(df$V1 == chr & df$V4 >= start & df$V4 <= end)
write.table(df[i,2], file=out, row.names=F, col.names=F, quote=F)
