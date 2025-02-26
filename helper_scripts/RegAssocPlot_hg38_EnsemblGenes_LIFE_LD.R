### creates RA plots based on hg38 as reference
  
createInputfilesRegAssocPlot1kg38Ensembl <- function(snpsubdata, r_on_server =T, leadsnp = "min_pval", resultdir = "auto",  margin = 500000, gene_validation_level = c("NOVEL",   "KNOWN"  ,  "PUTATIVE"), r2_objekt = "auto", path_ldreference = paste0(basicpath, "/01_daten/2306_lifea1_LDref_TOPMed/results/"), path_gen_annotation = "/net/ifs1/san_projekte/projekte/genstat/07_programme/rtools/1404_regional_assoc_plot_WO_broad/genes_hg38_ensembl/", file_gen_annotation = "Homo_sapiens.GRCh38.109_known_genes.RData", path_recombination = "/net/ifs1/san_projekte/projekte/genstat/01_daten/2306_recomb_rates_hg38_download/recomb-hg38/", file_recombination = "genetic_map_GRCh38_merged.tab", path_rprofile = paste0(path_regplot, "/../RProfile_hk/Rprofile_hk_150122hk.R"))
{
 
  ##### minidocumentation
  
  ## resultdir =  'auto' creates  a folder and subfolder named using the leadsnp inthe current working director to store created files required for the plot
  ## leadsnp = 'auto' usese the min p value snps of the file, otherwise specify a snp
  ## margin for some reason, BROAD imported in the gene file and the rrecombination file a range + - 500 000 bp . So did I 
  ## doplink calculate LD by plink, when the same data is analysed a sceond time FALSE helps to save time
  ## gene_validation_level is according to UCSC Table hg18 refSeqStatus, see below


  #   ############################################
  ### INIT
  
  snpsubdata = data.frame(snpsubdata)
  
  #### check input data
  
  stopifnot(all(c("snp", "position", "pvalue", "chr") %in% names(snpsubdata) ))
  snpsubdata = snpsubdata[ ,c("snp", "chr", "position", "pvalue")]
  
  if(any(grepl("X", snpsubdata$chr, ignore.case = T))) stop("do not use chr X as code for sex chr but chr 23")
  
  start = min(snpsubdata$position)
  stopifnot(is.numeric(start) & length(start) == 1)
  ende = max(snpsubdata$position)
  stopifnot(is.numeric(ende) & length(ende) == 1)
  chr = unique(as.numeric(snpsubdata$chr))
  stopifnot(is.numeric(chr) & length(chr) == 1 &(is.na(chr)==F))
  checke_NAs = sum(sapply(snpsubdata, is.na))
  if(checke_NAs !=0) stop("currently no NAs allowed in data.frame")
  
  ## Leadsnp
  
  if(leadsnp == "min_pval") leadsnp = snpsubdata[ order(snpsubdata$pvalue), "snp"][1] 
  stopifnot(leadsnp %in% snpsubdata$snp)
  
  
  ## directories 
  if(r_on_server ==T) {
    basicpath = "/net/ifs1/san_projekte/projekte/genstat/"
  } else basicpath = "R:/genstat/"
  
  
  path_regplot = paste0(basicpath, "07_programme/rtools/1404_regional_assoc_plot_WO_broad")
  
  directoryname = paste0("_files4RegionAssocPlot/_files4RegionAssocPlot_", leadsnp)
  directoryname = str_replace_all(directoryname, ":|-|>|<|\\.", "_")
  pos_leadsnp = snpsubdata[ snpsubdata$snp == leadsnp,"position"]
  
  if(identical(r2_objekt , "auto")) {    
    
    ##################################################
    ## LD calculation
    ### build plink command
    message("\n--------------------------\n")
    message("\nImporting LD Reference...\n")
    
    ### checke presence of leadsnp in 1kg
    ## snp annotaion hapmap
    require(data.table)
    
    ld = fread(paste0(path_ldreference, "LIFE-Adult_TOPMed_HG38_LD_chr",chr,".vcor.gz"), nThread = 5)
    checkin = ld[(ID_A == leadsnp) | (ID_B == leadsnp)]
    if (nrow(checkin) >0) checkin = T else checkin = F 
    message(" lead snp in Reference:", checkin, "\n")
    
    
    if(checkin) {
      ## LD 
      ld_match = ld[(ID_A == leadsnp) & (ID_B %in% snpsubdata$snp) |
                      (ID_B == leadsnp) & (ID_A %in% snpsubdata$snp)]
      
      if (nrow(ld_match) > 0){
        for(i in 1:nrow(ld_match)){
          if (ld_match[i,"ID_A"] == leadsnp) ld_match[i,"ID_A"] = ld_match[i,"ID_B"]
        }
      }
      
      ## file building
      leadsnp %in% ld$ID_A   #might not always be the right way to treat the leadsenp label
      
      data_tmp = data.frame(NAME = snpsubdata$snp,  CHR =  chr, PVAL = snpsubdata$pvalue, POS = snpsubdata$position)
      
      data_tmp = merge(data_tmp, ld_match[,c('ID_A', 'PHASED_R2' )], by.x="NAME", by.y = "ID_A", all.x = T, sort = F)
      data_tmp = data_tmp[, c("NAME", "CHR", "PVAL", "POS", "PHASED_R2")]
      
      #check if SNPs with no RSQR have an entry in the reference
      #if so set their LD to 0 instead of NA
      snps_no_ld = data_tmp["NAME"][is.na(data_tmp["PHASED_R2"]),]
      snps_in_reference = unique(c(ld$ID_A, ld$ID_B))
      snps_no_ld_in_reference = intersect(snps_no_ld, snps_in_reference)
      data_tmp["PHASED_R2"][data_tmp$NAME %in% snps_no_ld_in_reference,] = 0
      #ht(data_tmp,1)
      
      #data_tmp = rename(data_tmp, c(R2 = "RSQR"))
      names(data_tmp)[which(colnames(data_tmp)=="PHASED_R2")] = "RSQR"
    } else {
      message("As ID of lead SNP is not in Reference, no LD could be calculated...")
      data_tmp = data.frame(NAME = snpsubdata$snp,  CHR =  chr, PVAL = snpsubdata$pvalue, POS = snpsubdata$position, RSQR = NA)
    }
  } else data_tmp = r2_objekt
  ## save
  data_tmp_fn = paste0(directoryname, "/data_tmp")
  
  
  #########################################################
  ##  Recombination rate 
  
  message("\n--------------------------\n")
  message("\nImporting Recombination Rates ...\n")
  
  recomb_fn = paste0(path_recombination, file_recombination)
  recomb_fn
  
  recomb = fread(recomb_fn)
  
  setnames(recomb, names(recomb), c("chromosome","position", "comb_rate_cMPerMb", "map_cM"))
  recomb[,chromosome := str_replace(chromosome, "chr", "")]
  rate = data.frame(recomb[(chromosome == chr) & (position < ende + margin) & (position > start -  margin),  c("position", "comb_rate_cMPerMb"), with = F])
  ht(rate,2)  

  names(rate)[which(colnames(rate)=="position")] = "POS"
  names(rate)[which(colnames(rate)=="comb_rate_cMPerMb")] = "THETA"
  rate$DIST = abs(rate$POS -  pos_leadsnp)
  ## save
  rate_fn = paste0(directoryname, "/rate_tmp")
  
  ht(rate)
  
  #############################################################
  ##  Gene file
  message("\n--------------------------\n")
  message("\nImporting Genes...\n")
  
  require(biomaRt)
  
  
  loaded33 = load(paste0(path_gen_annotation, file_gen_annotation))
  loaded33
  require(data.table)
  
  genes2[,sort(unique(CHROMOSOME))]
  
  genes2[ CHROMOSOME == "X", CHROMOSOME :="23"]

  genes2 = genes2[ VALIDATION %in% gene_validation_level ]
  
  genes3 = genes2[ CHROMOSOME == chr & ((START > start - margin & START < ende +  margin) | 
                                          (STOP > start - margin & STOP < ende +  margin) |
                                          (START <= start - margin & STOP >= ende +  margin)
  )]
  
  
  setorder(genes3, STRAND, GENE, START, STOP)
  genes3[, start_nachhergen := c(START[-1], NA), by = GENE]
  genes3[, next_is_same_gene := ifelse(start_nachhergen <= STOP, T, F)]
  genes3[, end_vorhergen := c(NA, STOP[-length(STOP)]), by = GENE]
  genes3[, previous_is_same_gene := ifelse(end_vorhergen >= START, T, F)]
  
  genes3[, same_gene := next_is_same_gene | previous_is_same_gene]
  genes3[, rownum := 1:nrow(genes3)]

  genes3[, pseudo_id := ifelse(same_gene ==F | is.na(same_gene), rownum, GENE)]
  genes3 = genes3[, .(GENE ,CHROMOSOME, STRAND,    START    , STOP,VALIDATION    ,    biotype , pseudo_id)]
  setDF(genes3)

  genes4 = ddply(genes3, .variables=c("pseudo_id", 'STRAND'), .fun= function(df){
    
    df$START = min(df$START)
    df$STOP = max(df$STOP)
    df$VALIDATION = paste(unique(sort(df$VALIDATION)), collapse = ", ")
    unique(df)
    
  })
  
  if(any(genes4$STRAND) %nin% c(1,-1)) warning(paste0("found genes strand not only 1 and -1 as expected in ensembl..:", paste(unique(genes4$STRAND), collapse = " ")))
  genes4$STRAND = ifelse(genes4$STRAND == 1, "+", "-")
  
  ## save
  gene_fn = paste0(directoryname, "/genes_tmp")

  genes4
  message("\n----------------------\nDONE....looking forward to see the plot.....\n")
  
  input = c()
  input$myLocus = data_tmp[,c( "CHR","POS","NAME", "PVAL",  "RSQR")]
  input$myMap = rate
  input$myGenes =  genes4
  input$leadsnp = leadsnp
  input
}



customASplot_woBROAD = function (locus, map, genes, lead_snp = "auto",shownregion_kb = 500,flanking = 1000, best.pval = NULL, sf = c(4, 5), logpmax =  "auto", pch = 21, subtitle="", maintitle = "Regional Association Plot", gene_lines = "auto",dolegend = T, title_size = 1.1, subtitle_size = 0.8, print_lead_snp_name = F,print_lead_snp_pval =T, legendpos = "topleft", yticks = "auto", ylabel = "auto", lead_snp_name = "auto", textverschieb_standard = 5.5, center_lead_snp = T,  strongR2 = 0.8, moderateR2=0.5, weakR2 = 0.2 , plot_recrate_ontop =F, hard_xlim = F, col_recombirate = "lightblue", prefix_pWert ="P=",  usePVALcolumn_withoutlog10 = F, cex_genname=0.7, lead_snp_label_pos = 4, col_topSNP = "blue", ...) {
  
  
  if (dim(locus)[2] != 5 || is.null(dim(locus))) stop("Error: locus object has wrong dimension.")
  if (dim(map)[2] != 3 || is.null(dim(map))) stop("Error: map object has wrong dimension.")
  
  #1a. check column names 
  if (sum(c("CHR","POS","NAME","PVAL","RSQR") %in% colnames(locus)) != 5) 
    stop("Error: wrong column name in locus object.")
  if (sum(c("POS","THETA","DIST") %in% colnames(map)) != 3) 
    stop("Error: wrong column name in map object.")
  if (sum(c("START","STOP","STRAND","GENE", "VALIDATION", "CHROMOSOME") %in% colnames(genes)) != 6) 
    stop("Error: wrong column name in genes object.")
  
  # transform data if asked for ----
  if(usePVALcolumn_withoutlog10==F) locus$PVAL_transformed = -log10(locus$PVAL) else locus$PVAL_transformed = locus$PVAL 
  if(identical(logpmax,"auto")==T) logpmax =  max(locus$PVAL_transformed)
  
  
  #1b. check if given logpmax is suitable for given data----
  if (logpmax < max(locus$PVAL_transformed))
    warning("logpmax = ", logpmax,": largest value used for y axis  ",max(locus$PVAL_transformed))
  
  #2. define some needed variables
  ## identify best SNP and given chromosome
  if(lead_snp == "auto") bestSnp<-locus$NAME[locus$PVAL_transformed == max(locus$PVAL_transformed)] else  bestSnp<- lead_snp
  bestSnpDetails <- locus[locus$NAME == bestSnp,]
  
  ### limit to window of interest
  bestSnppos = bestSnpDetails$POS
  locus = locus[ (locus$POS < bestSnppos + shownregion_kb*1000) & (locus$POS > bestSnppos - shownregion_kb*1000),]
  ###
  
  
  actChr <- locus$CHR[1]
  ## identify best p-value if not given as parameter
  if (is.null(best.pval)) 
    best.pval <- bestSnpDetails$PVAL_transformed
  
  ## calculate plot range x axis
  # center_lead_snp = T: lead SNP is in the center, even if one side is then almost empty
  # center_lead_snp = F: first and last SNP in locus are used to set borders, lead SNP not necessarily in the center
  if(center_lead_snp == T){
    center <- bestSnppos
    lBorder <- center - shownregion_kb*1000 - flanking
    rBorder <- center + shownregion_kb*1000 + flanking
    xLength <- rBorder - lBorder
  }
  if(center_lead_snp == F){
    lBorder <- min(locus$POS) - flanking
    rBorder <- max(locus$POS) + flanking
    xLength <- rBorder - lBorder
    center <- lBorder + xLength/2
  }
  
  
  ## calculate plot range y axis
  # offset: space for gene names at bottom of plot 
  offset <- logpmax/sf[1]
  ylimPValue <- logpmax + offset
  # /4 is scaling the y axis of Recombination rate
  ylimRecomb <- ylimPValue/2
  
  #yadj is scaling difference between the two y axes
  yadj <- -offset + ylimPValue/sf[2]
  
  #3. shrink given data sets to plot boundaries
  keepMap <- subset(map, POS > lBorder & POS < rBorder)
  genes <- subset(genes, (START > lBorder & START < rBorder) | (STOP > lBorder & STOP < rBorder) | (START <= lBorder & STOP >= rBorder))
  
  #4. start plotting the image
  ## set margins of plot
  par(mar = c(3, 3.5, 2.5, 3.5)) # TODO param 5
  par(cex.axis = 0.75)
  ## calculate points for recombination rates and plot it as light blue line
  ## x coordinates as given in map object
  ## y coordinates have to be adjusted
  summary(keepMap$THETA/ylimRecomb)
  maxTheta  = max(keepMap$THETA)
  
  xy <- xy.coords(keepMap$POS, yadj + keepMap$THETA/(maxTheta/ylimRecomb))
  if(length(as.numeric(hard_xlim)) == 2)  plot(xy$x, xy$y, type = "l", col = col_recombirate, lwd = 1, ylim = c(-offset, logpmax), xlab = "", ylab = "", axes = F, xlim = hard_xlim) else plot(xy$x, xy$y, type = "l", col = col_recombirate, lwd = 1, ylim = c(-offset, logpmax), xlab = "", ylab = "", axes = F)
  
  ## simple box around plot
  box()
  ## labels for x axis (I want to use min/max of given region, without flanking)
  p1 <-  if(length(as.numeric(hard_xlim)) == 2) hard_xlim[1] else lBorder + flanking
  p5 <- if(length(as.numeric(hard_xlim)) == 2) hard_xlim[2] else rBorder - flanking
  p3 <- (p1 + p5)/2
  p2 <- (p1 + p3)/2
  p4 <- (p3 + p5)/2
  posXLabels <- c(p1, p2, p3, p4, p5)
  ## plot all three axes
  ## x axis
  axis(1, at = posXLabels, labels = format(round(posXLabels),big.mark=",",scientific=F), las = 1)
  mtext(paste("Chromosome", actChr, "(bp)", sep = " "), side = 1, line = 2, cex = 0.75)
  ## y axis left observed p
  ticks = if(identical(yticks, "auto")) signif(logpmax / 5,2) else {
    if(yticks > logpmax)
      warning("yticks is ", yticks, " that is larger as largest value in PVAL column (after transformation):",logpmax, " , consider lower values for yticks (i.e. the distance between labels on the y axis)")
    yticks}
  axis(2, at = seq(0, logpmax, ticks), labels = seq(0, logpmax, ticks), las = 1)
  
  ylabel = if(identical(ylabel, "auto")) "-log10(Observed p)" else ylabel
  
  breiteBeschriftungszahlen_yaxe = stringr::str_length(ticks)
  breiteZahlen_da = breiteBeschriftungszahlen_yaxe >2
  if(breiteZahlen_da==T) mtext(ylabel, side = 2, at = logpmax/2, line = breiteBeschriftungszahlen_yaxe-0.5, cex = 0.75) else mtext(ylabel, side = 2, at = logpmax/2, line = 2.5, cex = 0.7)
  ## y axis right recombination rate 
  ticks_recomb = ceiling(max(keepMap$THETA) / 5)
  
  
  
  axis(4, at = yadj + seq(from=0,to=100,by=ticks_recomb) / (maxTheta/ylimRecomb), labels = seq(from=0,to=100,by=ticks_recomb), las = 1, col.axis = "blue")
  mtext("Recombination rate (cM/Mb)", side = 4,  line = 2, col = "blue", cex = 0.75)
  
  ## mark 0 of left y axis as dotted line
  lines(c(lBorder, rBorder), c(0, 0), lty = "dotted", lwd = 1, col = "black")
  
  ## Attention: best snp is plotted at position logpmax, if range of y axis is not large enough
  if (unique(best.pval) > logpmax) posBestSnp <- logpmax  else posBestSnp <- best.pval
  points(bestSnpDetails$POS, posBestSnp, pch=pch, cex=2.5, bg=col_topSNP)
  
  if(identical(lead_snp_name, "auto")) mylead_snp_name = bestSnpDetails$NAME else mylead_snp_name = lead_snp_name
  if(print_lead_snp_name) text(bestSnpDetails$POS, posBestSnp, labels = mylead_snp_name, pos=lead_snp_label_pos, offset=1)
  if(print_lead_snp_pval) text(bestSnpDetails$POS, posBestSnp, labels = c(paste(prefix_pWert,signif(bestSnpDetails$PVAL,3), sep="")), pos=lead_snp_label_pos, offset=1)
  
  ## divide locus object in different subsets depending on r^2
  strong.ld <- subset(locus, locus$NAME != bestSnp & locus$RSQR >= strongR2)
  moderate.ld <- subset(locus, locus$NAME != bestSnp & locus$RSQR >= moderateR2 & locus$RSQR < strongR2)
  weak.ld <- subset(locus, locus$NAME != bestSnp & locus$RSQR >= weakR2 & locus$RSQR < moderateR2)
  not.in.ld <- subset(locus, locus$NAME != bestSnp & locus$RSQR < weakR2)
  na.ld <- subset(locus, locus$NAME != bestSnp & is.na(locus$RSQR))
  
  ## define colors and plot r^2 values
  colors <- c("red", "orange", "yellow", "grey", "white")
  
  ## plot legend of r^2
  ltext <- rbind(paste0("[",strongR2,"-1.0]"), paste0("[",moderateR2,"-",strongR2,")"), paste0("[",weakR2,"-",moderateR2,")"), paste0("[0.0-",weakR2,")"), "NA")
  
  if(dolegend) legend(legendpos, legend = ltext, title = expression(r^2), fill = colors, cex = 0.7, ...)
  points(na.ld$POS, na.ld$PVAL_transformed, pch = pch, cex = 1, bg = colors[5])
  points(not.in.ld$POS, not.in.ld$PVAL_transformed, pch = pch, cex = 1, bg = colors[4])
  
  points(weak.ld$POS, weak.ld$PVAL_transformed, pch = pch, cex = 1.1, bg = colors[3])
  
  points(moderate.ld$POS, moderate.ld$PVAL_transformed, pch = pch, cex = 1.1, bg = colors[2])
  
  points(strong.ld$POS, strong.ld$PVAL_transformed, pch = pch, cex = 1.1, bg = colors[1])
  
  ## plot of genes with name and strand orientation within the study region
  ## how much space is there for writing gene names?
  geneLines <- if(gene_lines == "auto") ceiling(offset/(ticks/2)) else gene_lines
  message(geneLines,fill=T)  
  ## if there are more genes than gene lines in the plot, order genes ascending by start position of gene
  ## not perfect, but works quite good to avoid into each other plotted genes/gene names
  if (nrow(genes) > geneLines) {
    genes<-genes[order(genes$START),]
  }
  if(dim(genes)[1]>0) {
    for (i in 1:nrow(genes)) {
      actGene <- genes[i, ]
      geneName <- as.character(actGene$GENE)
      #message("-",geneName,fill=T)
      geneStart <- actGene$START
      geneStop <- actGene$STOP
      center <- (geneStart + geneStop)/2
      if(center <  lBorder)  center  = geneStop
      if(center >  rBorder)  center = geneStart
      
      adj <- -offset + (ticks/2)*(i%%geneLines) -ticks/4.5
      
      ## print arrow for gene, depending on strand information
      strandInfo <- actGene$STRAND
      if (strandInfo == "+") 
        arrows(geneStart, adj, geneStop, adj, length = 0.05, lwd = 2, code = 2, lty = "solid", col = "darkgreen")
      else arrows(geneStart, adj, geneStop, adj, length = 0.05, lwd = 2, code = 1, lty = "solid", col = "darkgreen")
      
      ## print gene names
      textverschieb = ticks/textverschieb_standard
      if (!is.na(geneName)) text(center, adj + textverschieb, labels = geneName, cex = cex_genname)
      
    }
  }
  ## plot title and subtitle
  title(main=maintitle, cex.main = title_size)
  mtext(subtitle, cex = subtitle_size)
  if(plot_recrate_ontop ==T)   points(xy$x, xy$y, type = "l", col = col_recombirate, lwd = 1, ylim = c(-offset, logpmax), xlab = "", ylab = "", axes = F)
}
