liftOVerHg38TOHg19 <- function (tolift, try_unlifted_via_ensembl = T, return_as_data_table = F, path_chainfile = paste0(basicpath, "/07_programme/rtools/150303_liftover_via_R/data/hg38ToHg19.over.chain")) {
  message('Using function liftOVerHg38TOHg19()...')
  ### input: tolift : data.table or data.frame with hg38 chr, pos und ID
  ###        try_unlifted_via_ensembl: Check unlifted SNPs with Ensembl
  ###         return_as_data_table: return data.table or data_frame
  ###         path_chainfile: path to chainfile
  
  ##packages
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(reshape2)
  library(stringr)
  if (try_unlifted_via_ensembl ) library(biomaRt)
  ## checken
  stopifnot(c('snps', 'chr', 'pos') %in% names(tolift))
  
  tolift_ori = tolift
  
  ## data.tablen
  tolift = data.table(tolift)
  check1 = tolift[,stopifnot(all(na.omit(chr) %in% 1:25))]
  check1
  tolift[,chr := as.character(chr)]
  
  
  ## filter NAs
  nafilter = is.na(tolift$chr) | is.na(tolift$pos)
  # mytable(nafilter)
  message('Not lifting ', length(unique(tolift[ nafilter, snps])), " SNPs since NA was given in chr or pos...")
  tolift_na = tolift[ nafilter]
  tolift = tolift[ nafilter ==F]
  
  
  ## convert to genomic range
  
  ### chr umkodieren
  if(any(tolift$chr %in% as.character(23:25))) message("I assume chromosome 23 == X, 24 == Y, and 25 == M")
  tolift[ chr =='23', chr := "X"]
  tolift[ chr =='24', chr := "Y"]
  tolift[ chr =='25', chr:= "M"]
  
  tolift[ , chr := paste0("chr",chr)]
  
  
  tolift_ir  = IRanges(start = tolift$pos, end = tolift$pos)
  tolift_ir
  
  tolift_gr = GRanges(seqnames = Rle(tolift$chr), 
                      ranges = IRanges(start = tolift$pos, end = tolift$pos),
                      strand = Rle("*", length(seqnames))
  )
  mcols(tolift_gr)$snps = tolift$snps   
  tolift_gr    
  
  
  
  ## ----dolift--------------------------------------------------------------
  ch = import.chain(path_chainfile)
  ch
  names(ch)
  lifted = liftOver(tolift_gr, ch)
  class(lifted)
  lifted
  lifted = unlist(lifted)
  lifted
  
  genome(lifted) = "hg19"
  
  names(lifted@elementMetadata) = "snps"
  ## convert to data table and merge
  lifted_dt = data.table(snps = lifted$snps, chr_hg19 = as.character(seqnames(lifted)), pos_hg19 = ranges(lifted)@start)
  
  setkey(lifted_dt, snps)
  setkey(tolift, snps)
  
  res = lifted_dt[tolift]
  # res
  
  res = res[ ,chr_hg19:= str_replace(chr_hg19, "chr", "")]
  res = res[chr_hg19== "X", chr_hg19 := '23']
  res = res[chr_hg19== "Y", chr_hg19 := '24']
  res = res[chr_hg19== "M", chr_hg19 := '25']
  
  res = res[ ,chr:= str_replace(chr, "chr", "")]
  res = res[chr== "X", chr := '23']
  res = res[chr== "Y", chr := '24']
  res = res[chr== "M", chr := '25']
  res[,liftmethod := "liftover"]
  # res
  
  
  
  ## ----vgl--------------------------------------------------------------
  lifted_sucessful = res[is.na(pos_hg19)==F]
  
  lifted_failed = res[is.na(pos_hg19) ]
  lifted_failed
  
  ## redo failed liftings via USC
  message("Could not lift via 'hg38ToHg19.over.chain' ", dim(lifted_failed)[1], "  SNPs...")
  
  if(try_unlifted_via_ensembl & nrow(lifted_failed) >0) {
    
    
    message("Try getting lifting-failed positions online via Ensemble...")

    
    
    ensemblsnp =  try(useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp",host = "feb2014.archive.ensembl.org", path = "/biomart/martservice", archive = FALSE))
    
 
    ensemblsnp
      
      mycategs = c("refsnp_id", "chr_name", "chrom_start" )
      snpattrib =  getBM(attributes=mycategs, filters = 'snp_filter', values = lifted_failed$snps, mart = ensemblsnp)
      
      message("got information for ", length(unique(snpattrib$refsnp_id)), " SNPs...")
      
      snpattrib = data.table(snpattrib)
      snpattrib[,chrlaenge := str_length(chr_name)]
      kurz = snpattrib[ chrlaenge <3]
      lang = snpattrib[ chrlaenge >2]
      
      save2kill = lang[ refsnp_id %in% kurz$refsnp_id, refsnp_id]
      lang = lang[ refsnp_id %in% save2kill ==F]
      lang
      
      snpattrib = rbindlist(list(kurz, lang))
      snpattrib[ ,chrlaenge :=NULL]
      
      setnames(snpattrib,c("refsnp_id", "chr_name", "chrom_start"), c("snps", "chr_hg19", "pos_hg19"))
      snpattrib[ chr_hg19 == "X", chr_hg19 :='23']
      snpattrib[ chr_hg19 == "Y", chr_hg19 :="24"]
      snpattrib[ chr_hg19 == "MT", chr_hg19 :="25"]
      
      
      
      snpattrib[ , liftmethod := "pos_via_ID_ensemblehg19"]
      
      tomerge = lifted_failed[,list(snps, chr, pos)]
      setkey(tomerge, snps)
      setkey(snpattrib, snps)
      
      snpattrib2 = tomerge[ snpattrib]
      setcolorder(snpattrib2, names(lifted_sucessful))
      lifted_sucessful = rbindlist(list(lifted_sucessful, snpattrib2))
      lifted_failed = lifted_failed[ snps %in% lifted_sucessful$snps ==F]
      
      res = rbindlist(list(lifted_sucessful, lifted_failed))
      
  #}
  }
  
  ## add NA SNPs
  res
  tolift_na[,chr_hg19:=NA]
  tolift_na[,pos_hg19:=NA]
  tolift_na[,liftmethod:="notlifted_(chr_and/or_pos_was_NA)"]
  setcolorder(tolift_na, names(res))
  
  res  = rbindlist(list(res, tolift_na))
  message("Succesfully lifted: " ,dim(lifted_sucessful)[1], " entries...\n")
  message("Lifting failes: " ,dim(lifted_failed)[1], " entries...\n")
  message("Not lifted (pos and/or chr was missing): " ,dim(tolift_na)[1], " entries...\n")
  
  message("\n returning following chromosomes: ", paste(sort(unique(res$chr_hg19)), collapse = " "))
  
  stopifnot(all(tolift_ori$snps %in% res$snps))
  message('FINISHED function liftOVerHg38TOHg19()...')
  
  if(return_as_data_table) return(res) else return(data.frame(res))
}
