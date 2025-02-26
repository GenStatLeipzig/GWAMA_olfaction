suppressPackageStartupMessages(library(data.table, quietly=T, warn.conflicts=F, verbose=F))
suppressPackageStartupMessages(library(coloc, quietly=T, warn.conflicts=F, verbose=F))

args <- commandArgs(trailingOnly=TRUE)

region <- args[1]
qtl.file <- args[2]
chr <- args[3]
N.qtl <- args[4]
out.file <- args[5]
cs.prefix <- args[6]
label <- args[7]

# load olfaction and qtl data
olf <- as.data.frame(fread("coloc_data_cropped_regions_qc_filtered.csv.gz", header=T))
gwas <- olf[which(olf$region == region),c(1:8,39,42:44,48,49)]
qtl <- as.data.frame(fread(qtl.file, header=T))

# lift qtl dat
lift <- read.delim(paste("lift/qtl.mapping.chr",chr,".vcfliftover.txt", sep=""), header=F)
names(lift) <- c("chr", "pos38", "pos37")
m <- match(qtl$BP, lift$pos37)
qtl$BP <- lift$pos38[m]

qtl$N <- N.qtl

# in both gwas and qtl, the NEA / EA allele is random with respect to ref / alt
# in qtl, dont use A1, it's * and won't match to gwas
### in qlt,  A1==the effect (coded) allele, A2==the other allele
### in gwas, ea==the effect (coded) allele, aa==the other allele

qtl <- cbind(qtl, variant_id=paste(qtl$Chr, qtl$BP, qtl$A1,qtl$A2, sep="_"))
gwas <- cbind(gwas, variant_id=paste(gwas$chrom, gwas$pos, gwas$ea, gwas$aa, sep="_"))
merged_1 <- merge(gwas, qtl, by="variant_id", all=FALSE)

# try flipped alleles in qtl
qtl$variant_id <- paste(qtl$Chr, qtl$BP, qtl$A2,qtl$A1, sep="_")
merged_2 <- merge(gwas, qtl, by="variant_id", all=FALSE)
cat("is there any overlap between the two merged dfs?", any(unique(merged_1$BP) %in% unique(merged_2$BP)), "\n")


# in merged_2, we need to flip the beta/p of qtl so that we measure the effect for the same alleles
merged_2$b <- -(merged_2$b)

# combine both merged files and reorder
merged <- rbind(merged_1, merged_2)
o <- order(merged$pos)
merged <- merged[o,]


# how many metal variants could not be matched to gtex?
incl <- length(which(gwas$pos %in% unique(merged$pos)))
total <- length(unique(gwas$pos))
cat(incl, "out of", total, "metal variants, (", round((incl/total)*100), "%) could be matched to gtex variants\n")

i <- which(merged$Freq > 0.5)
if (length(i)>0){
    merged$Freq[i] <- 1-merged$Freq[i]
}

resall<-NULL
for (phenotype in unique(merged$phenotype)){
    cat("running phenotype", phenotype, "\n")
    for (gene in unique(merged$Probe)) {
        cdat <- merged[which(merged$Probe==gene & merged$phenotype==phenotype),]
        # has to overlap either start or end or lie in the region

        # there were some cases, where the lifting assigned one b37 to two b38 positions, creating duplicate ids.
        # they throw an error in the coloc remove them but then check manually how many there
        if (length(unique(cdat$markerID)) < nrow(cdat)){
            tab <- as.data.frame(table(cdat$markerID))
            ids <- tab$Var1[which(tab$Freq > 1)]
            # remove whole snp because I don't know which to keep
            cat("CustomWarning: removed SNPs", ids, "because they were duplicates in the df\n")
            for (id in ids){
                print(cdat[which(cdat$markerID==id), ])
                cdat <- cdat[-which(cdat$markerID==id), ]
            }
        }
        if (nrow(cdat) >= 30){
            # according to Claudia the threshould should be 50, but put a bit lower here and filter later
            # make lists of both GWAS sets for coloc
            #D1=gwas
            D1 <- list(type = "quant", beta=cdat$betaFEM, varbeta=cdat$seFEM^2, pvalues=cdat$pFEM, N=cdat$totalN, MAF=cdat$nWeightedMAF, snp=cdat$markerID, sdY=1)
            #D2=qtl
            D2 <- list(type = "quant", beta=cdat$b, varbeta=cdat$SE^2, pvalues=cdat$p, N=cdat$N, MAF=cdat$Freq, snp=cdat$markerID, sdY=1)
            coloc.res <- coloc.abf(D1, D2)
            res <- as.data.frame(t(as.data.frame(coloc.res$summary)))
            res$probe <- gene
            res$region <- region
            res$phenotype <- phenotype
            res$qtl <- label
            res$N.qtl <- N.qtl
            # add the most likely causal variant
            df <- coloc.res$results
            i <- which(df$SNP.PP.H4 ==  max(df$SNP.PP.H4))[1]
            res$max.SNP.PP.H4 <- df$SNP.PP.H4[i]
            res$ID.SNP.PP.H4 <- df$snp[i]
            resall <- rbind(resall,res)
            # only if coloc is sign.
            if(res$PP.H4.abf > 0.8){
                # save 95% conf set
                o <- order(df$SNP.PP.H4,decreasing=TRUE)
                cs <- cumsum(df$SNP.PP.H4[o])
                w <- which(cs > 0.95)[1]
                write.table(df[o,][1:w,]$snp, file=paste(cs.prefix,"_", gene,"_",phenotype,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
            }
        }
    }
}
if(!is.null(resall)){
    resall <- resall[order(resall$PP.H4.abf,decreasing=T),]
    write.table(resall, file=out.file, row.names=F, col.names=T,sep="\t", quote=F)
}else{
    cat("couldnt do any coloc tests, no data available\n")
}

