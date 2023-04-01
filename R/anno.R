#' @title calculate genomic control inflation factor for a QTL/GWAS summary statistics dataset.
#'
#' @param summaryDT A data.frame containing one or two columns: p-value (required) and group (optional)
#' @import data.table
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "https://raw.githubusercontent.com/dingruofan/exampleData/master/eqtl/MMP7_qtlDF.txt"
#' qtl <- data.table::fread(url1, sep="\t")
#'
#' # calculate lambda value with all variants
#' xQTLanno_calLambda(qtl[,.(pValue)])
#'
#' # calculate lambda value for each group:
#' qtl$groups <- sample(c(0,1),size = nrow(qtl), replace = TRUE)
#' xQTLanno_calLambda(qtl[,.(pValue, groups)])
#' }
xQTLanno_calLambda <- function(summaryDT){
  .<-NULL
  message("== Number of variants: ", nrow(summaryDT), "  ",date())
  if(ncol(summaryDT)==1){
    summaryDT <- summaryDT[,1]
    names(summaryDT) <- c("pval")
    summaryDT <- na.omit(summaryDT)
    data.table::setDT(summaryDT)
    message("calculating lamdba: ", date())
    # http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
    lamdba_value <- median(qchisq(1-summaryDT$pval,1))/qchisq(0.5,1)
    message("== Lamdba: ", lamdba_value)
    lamdba_value <- data.table(groups=NA, lambda=lamdba_value)
  }else if(ncol(summaryDT)==2){
    summaryDT <- summaryDT[,1:2]
    names(summaryDT) <- c("pval", "groups")
    summaryDT <- na.omit(summaryDT)
    data.table::setDT(summaryDT)
    lamdba_value <- summaryDT[,.(lambda= median(qchisq(1-.SD$pval,1))/qchisq(0.5,1)),by="groups"]
    message("== Lamdba: ", lamdba_value)
  }
  return(lamdba_value)
}


#' @title enrichment analysis for GWAS / QTL signals in functional elements, including enhancer, promoter, CPG, and TFs
#'
#' @param snpInfo A data.table/data.frame with two columns: chromosome and position.
#' @param genomeVersion "hg38" (default) or "hg19". Note: hg19 will be converted to hg38 automatically.
#' @param enrichElement "Promoter", "Enhancer" or "TF".
#' @param distLimit Defaults: 1e6
#' @importFrom GenomicRanges distanceToNearest mcols
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "https://github.com/dingruofan/exampleData/raw/master/gwas/gwasSub.txt.gz"
#' snpInfo <- data.table::fread(url1, sep="\t")
#' enrichHits <- xQTLanno_enrich(snpInfo,enrichElement="Enhancer")
#' }
xQTLanno_enrich <- function(snpInfo="",  genomeVersion="hg38", enrichElement="Promoter", distLimit=1e6){
  chrom <- pValue <- tx_name <- distance <- V1 <- V2 <- V3 <- NULL
  . <- NULL

  message("==> Start checking variants...")
  # check snpinfo:
  snpInfo <- data.table::as.data.table(snpInfo[,1:3])
  snpInfo <- na.omit(snpInfo)
  names(snpInfo) <- c("chrom", "pos", "pValue")
  snpInfo$chrom <- as.character(snpInfo$chrom)
  snpInfo$pos <- as.integer(snpInfo$pos)
  snpInfo$pValue <- as.numeric(snpInfo$pValue)
  snpInfo <- snpInfo[chrom %in% paste0("chr",c(1:22,"X","Y","M"))]
  snpInfo <- snpInfo[pValue>0,]
  if(nrow(snpInfo)<1){
    stop("Number of variants < 1...")
  }
  # create snpinfo Grange object:
  snpRanges <- GenomicRanges::GRanges(snpInfo$chrom,
                                      IRanges::IRanges(snpInfo$pos, snpInfo$pos),
                                      strand= "*",
                                      snpInfo[,"pValue"]
  )
  snpInfoNew <- data.table::data.table()
  snpRanges_hg38 <- GenomicRanges::GRanges()
  if(genomeVersion == "hg19"){
    message("== Converting vairants' coordinate to hg38... ")

    if(suppressMessages(!requireNamespace("rtracklayer"))){
      message("Package [rtracklayer] is not installed! please install [rtracklayer] with following: ")
      message("---------")
      message('\"if (!require("BiocManager", quietly = TRUE)); BiocManager::install("rtracklayer")\"')
      message("---------")
    }

    path = system.file(package="xQTLbiolinks", "extdata", "hg19ToHg38.over.chain")
    ch = rtracklayer::import.chain(path)
    snpRanges_hg38 <- unlist(rtracklayer::liftOver(snpRanges, ch))
    # retain the SNPs located in chromosome 1:22
    snpRanges_hg38 <- snpRanges_hg38[which(as.character(GenomeInfoDb::seqnames(snpRanges_hg38)) %in% paste0("chr", c(1:22,"X","Y", "M"))),]
    message("== ",length(snpRanges_hg38),"/",nrow(snpInfo)," (",round(length(snpRanges_hg38)/nrow(snpInfo)*100,2),"%)"," left.")
    snpInfoNew <- data.table::data.table(chrom = as.character(GenomeInfoDb::seqnames(snpRanges_hg38)),
                                         pos=as.data.table(IRanges::ranges(snpRanges_hg38))$start,
                                         pValue=snpRanges_hg38$pValue)
    snpRanges <- data.table::copy(snpRanges_hg38)
    snpInfo <- data.table::copy(snpInfoNew)
    rm(path, ch, snpInfoNew, snpRanges_hg38)
  }
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  if( enrichElement == "Promoter"){
    # extract promoters:
    message("      enrichment analysis of promoters...")
    promoterParts <- suppressWarnings(GenomicFeatures::promoters(txdb))
    promoterParts <- promoterParts[which(as.character(GenomeInfoDb::seqnames(promoterParts)) %in% paste0("chr", c(1:22,"X","Y", "M"))),]
    #
    # preDist <- as.data.table(GenomicRanges::distanceToNearest(snpRanges, promoterParts[subjectHits(precede(snpRanges, promoterParts, select="all")),]))
    # followDist <- as.data.table(GenomicRanges::distanceToNearest(snpRanges, promoterParts[subjectHits(follow(snpRanges, promoterParts, select="all")),]))
    # distPreFollow <- merge(preDist, followDist, by="queryHits", suffixes=c("_pre", "_follow"))
    nearestDist <- as.data.table(GenomicRanges::distanceToNearest(snpRanges, promoterParts))
    if(nrow(nearestDist)>0){
      nearestDist <- do.call(cbind, list(snpInfo[nearestDist$queryHits,],  as.data.table(GenomicRanges::mcols(promoterParts))[nearestDist$subjectHits,.(anno=tx_name)],nearestDist[,.(dist=distance)] ))
      nearestDist <- nearestDist[dist<distLimit,][order(dist)]
    }

  }else if(enrichElement == "Enhancer"){
    # import enhancer dataset:
    message("      enrichment analysis of enhancers...")
    enhancerAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "geneHancer.bed.gz"), header=FALSE)
    enhancerRanges <- GenomicRanges::GRanges(enhancerAnno$V1,
                                             IRanges::IRanges(enhancerAnno$V2, enhancerAnno$V3), strand= "*")
    nearestDist <- as.data.table(GenomicRanges::distanceToNearest(snpRanges, enhancerRanges))
    if(nrow(nearestDist)>0){
      nearestDist <- do.call(cbind, list(snpInfo[nearestDist$queryHits,], enhancerAnno[nearestDist$subjectHits,.(anno=paste0(V1,":",V2,"-",V3))], nearestDist[,.(dist=distance)] ))
      nearestDist <- nearestDist[dist<distLimit,][order(dist)]
    }

  }else if(enrichElement == "TF"){
    # import tf cluster dataset:
    message("      enrichment analysis of tf clusters...")
    tfAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "tf_filtered_sorted_merge_sub.bed.gz"), header=FALSE)
    tfRanges <- GenomicRanges::GRanges(tfAnno$V1,
                                       IRanges::IRanges(tfAnno$V2, tfAnno$V3), strand= "*")
    nearestDist <- as.data.table(GenomicRanges::distanceToNearest(snpRanges, tfRanges))
    if(nrow(nearestDist)>0){
      nearestDist <- do.call(cbind, list(snpInfo[nearestDist$queryHits,], tfAnno[nearestDist$subjectHits,.(anno=paste0(V1,":",V2,"-",V3))], nearestDist[,.(dist=distance)] ))
      nearestDist <- nearestDist[dist<distLimit,][order(dist)]
    }
  }
  return(nearestDist)
}


#' @title annotate all signals in GWAS / QTL dataset by genome location
#'
#' @param snpInfo A data.table/data.frame with three columns: chromosome, position, and P-value.
#' @param genomeVersion "hg38" (default) or "hg19". Note: hg19 will be converted to hg38 automatically.
#' @import data.table
#' @import stringr
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom GenomicFeatures intronicParts exonicParts transcripts
#' @importFrom IRanges findOverlaps
#'
#' @return A data.table object of variants' genomics distribution
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "https://github.com/dingruofan/exampleData/raw/master/gwas/gwasSub.txt.gz"
#' snpInfo <- data.table::fread(url1, sep="\t")
#' snpHits <- xQTLanno_genomic(snpInfo)
#' }
xQTLanno_genomic <- function(snpInfo="", genomeVersion="hg38"){
  chrom <- pValue <- V1 <- V2 <- V3 <- tx_name <- anno <- TXID <- tx_id <- type <- pos <- proportion <- NULL
  .<- NULL
  message("==> Start checking variants...")
  # check snpinfo:
  snpInfo <- data.table::as.data.table(snpInfo[,1:3])
  snpInfo <- na.omit(snpInfo)
  names(snpInfo) <- c("chrom", "pos", "pValue")
  snpInfo$chrom <- as.character(snpInfo$chrom)
  snpInfo$pos <- as.integer(snpInfo$pos)
  snpInfo$pValue <- as.numeric(snpInfo$pValue)
  snpInfo <- snpInfo[chrom %in% paste0("chr",c(1:22,"X","Y","M"))]
  snpInfo <- snpInfo[pValue>0,]
  if(nrow(snpInfo)<1){
    stop("Number of variants < 1...")
  }
  # create snpinfo Grange object:
  snpRanges <- GenomicRanges::GRanges(snpInfo$chrom,
                                      IRanges::IRanges(snpInfo$pos, snpInfo$pos),
                                      strand= "*",
                                      snpInfo[,"pValue"]
  )
  snpInfoNew <- data.table::data.table()
  snpRanges_hg38 <- GenomicRanges::GRanges()
  if(genomeVersion == "hg19"){
    message("== Converting vairants' coordinate to hg38... ")

    if(suppressMessages(!requireNamespace("rtracklayer"))){
      message("Package [rtracklayer] is not installed! please install [rtracklayer] with following: ")
      message("---------")
      message('\"if (!require("BiocManager", quietly = TRUE)); BiocManager::install("rtracklayer")\"')
      message("---------")
    }

    path = system.file(package="xQTLbiolinks", "extdata", "hg19ToHg38.over.chain")
    ch = rtracklayer::import.chain(path)
    snpRanges_hg38 <- unlist(rtracklayer::liftOver(snpRanges, ch))
    # retain the SNPs located in chromosome 1:22
    snpRanges_hg38 <- snpRanges_hg38[which(as.character(GenomeInfoDb::seqnames(snpRanges_hg38)) %in% paste0("chr", c(1:22,"X","Y", "M"))),]
    message("== ",length(snpRanges_hg38),"/",nrow(snpInfo)," (",round(length(snpRanges_hg38)/nrow(snpInfo)*100,2),"%)"," left.")
    snpInfoNew <- data.table::data.table(chrom = as.character(GenomeInfoDb::seqnames(snpRanges_hg38)),
                                         pos=as.data.table(IRanges::ranges(snpRanges_hg38))$start,
                                         pValue=snpRanges_hg38$pValue)
    snpRanges <- data.table::copy(snpRanges_hg38)
    snpInfo <- data.table::copy(snpInfoNew)
    rm(path, ch, snpInfoNew, snpRanges_hg38)
  }
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  ############## import datasets:
  message("==> Start importing annotations...")
  message("      cpg...")
  cpgAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "cpgLsland.bed.gz"), header=FALSE)
  cpgRanges <- GenomicRanges::GRanges(cpgAnno$V1,
                                      IRanges::IRanges(cpgAnno$V2, cpgAnno$V3), strand= "*")

  # import enhancer dataset:
  message("      enhancer...")
  enhancerAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "geneHancer.bed.gz"), header=FALSE)
  enhancerRanges <- GenomicRanges::GRanges(enhancerAnno$V1,
                                           IRanges::IRanges(enhancerAnno$V2, enhancerAnno$V3), strand= "*")

  # import tf cluster dataset:
  message("      tf cluster...")
  tfAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "tf_filtered_sorted_merge_sub.bed.gz"), header=FALSE)
  tfRanges <- GenomicRanges::GRanges(tfAnno$V1,
                                     IRanges::IRanges(tfAnno$V2, tfAnno$V3), strand= "*")

  # extract introns:
  message("      introns...")
  intronParts <-  GenomicFeatures::intronicParts(txdb)
  intronParts <- intronParts[which(as.character(GenomeInfoDb::seqnames(intronParts)) %in% paste0("chr", c(1:22,"X","Y", "M"))),]
  intronAnno <- as.data.table(intronParts)

  # extract exons:
  message("      exons...")
  exonParts <- GenomicFeatures::exonicParts(txdb)
  exonParts <- exonParts[which(as.character(GenomeInfoDb::seqnames(exonParts)) %in% paste0("chr", c(1:22,"X","Y", "M"))),]
  exonAnno <- as.data.table(exonParts)

  # extract promoters:
  message("      promoters...")
  promoterParts <- suppressWarnings(GenomicFeatures::promoters(txdb))
  promoterParts <- promoterParts[which(as.character(GenomeInfoDb::seqnames(promoterParts)) %in% paste0("chr", c(1:22,"X","Y", "M"))),]

  # all transcripts:
  txinfo <- as.data.table(transcripts(txdb))
  txinfo <- txinfo[which(seqnames %in% paste0("chr", c(1:22,"X","Y", "M"))),]

  ############ start find overlap with annotations:
  message("==> Start annotating variants...")
  cpgHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                              subject = cpgRanges,
                                                              maxgap = 1,
                                                              ignore.strand=TRUE ))
  cpgHits <- cbind(snpInfo[cpgHits$queryHits,],cpgAnno[cpgHits$subjectHits,.(anno=paste0(V1,":",V2,"-",V3), type="cpg")])
  message("      cpg hits: ",nrow(cpgHits))
  rm(cpgAnno, cpgRanges)

  # annotate snp with enhancers:
  enhancerHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                                   subject = enhancerRanges,
                                                                   maxgap = 1,
                                                                   ignore.strand=TRUE ))
  enhancerHits <- cbind(snpInfo[enhancerHits$queryHits,],enhancerAnno[enhancerHits$subjectHits,.(anno=paste0(V1,":",V2,"-",V3), type="enhancer")])
  message("      enhancer hits: ",nrow(enhancerHits))
  rm(enhancerAnno, enhancerRanges)

  # annotate snp with tf clusters:
  tfHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                             subject = tfRanges,
                                                             maxgap = 1,
                                                             ignore.strand=TRUE ))
  tfHits <- cbind(snpInfo[tfHits$queryHits,],tfAnno[tfHits$subjectHits,.(anno=paste0(V1,":",V2,"-",V3), type="tfCluster")])
  message("      tf hits: ",nrow(tfHits))
  rm(tfAnno, tfRanges)

  # annotate snp with introns:
  intronHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                                 subject = intronParts,
                                                                 maxgap = 1,
                                                                 ignore.strand=TRUE ))
  convertCharacter <- function(x){ a=as.list(x);paste0(unlist(a), collapse = ",") }
  if(nrow(intronHits)>0){
    intronHits <- cbind(snpInfo[intronHits$queryHits,],intronAnno[intronHits$subjectHits,.(anno=tx_name, type="exon")])
    intronHits$anno <- unlist(lapply(intronHits$anno, convertCharacter))
    # collapse annotation:
    intronHits <- intronHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  }else{
    intronHits <- cpgHits[0,]
  }
  message("      intron hits: ",nrow(intronHits))
  rm(intronParts, intronAnno)

  # annotate snp with exons:
  exonHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                               subject = exonParts,
                                                               maxgap = 1,
                                                               ignore.strand=TRUE ))
  if(nrow(exonHits)>0){
    exonHits <- cbind(snpInfo[exonHits$queryHits,],exonAnno[exonHits$subjectHits,.(anno=tx_name, type="exon")])
    exonHits$anno <- unlist(lapply(exonHits$anno, convertCharacter))
    # collapse annotation:
    exonHits <- exonHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  }else{
    exonHits <- cpgHits[0,]
  }
  message("      exon hits: ",nrow(exonHits))
  rm(exonParts, exonAnno)

  # annotate snp with promoters:
  promoterHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                                   subject = promoterParts,
                                                                   maxgap = 1,
                                                                   ignore.strand=TRUE ))
  promoterHits <- cbind(snpInfo[promoterHits$queryHits,], as.data.table(GenomicRanges::mcols(promoterParts))[promoterHits$subjectHits,.(anno=tx_name, type="promoter")])
  # collapse annotation:
  promoterHits <- promoterHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  message("      promoter hits: ",nrow(promoterHits))
  rm(promoterParts)

  # annotate snp with cds:
  cdsHits <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges,
                                                                                                txdb,
                                                                                                VariantAnnotation::CodingVariants())) ))
  if(nrow(cdsHits)>0){
    cdsHits <- cbind(snpInfo[cdsHits$QUERYID,], cdsHits[,.(tx_id=as.integer(TXID),type="cds")])
    cdsHits <- merge(cdsHits, txinfo[,.(tx_id,anno=tx_name)],by="tx_id")[,-c("tx_id")]
    # collapse annotation:
    cdsHits <- cdsHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  }else{
    cdsHits <- cpgHits[0,]
  }
  message("      cds hits: ",nrow(cdsHits))

  # annotate snp with 5 utr:
  utr5Hits <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges,
                                                                                                 txdb,
                                                                                                 VariantAnnotation::FiveUTRVariants())) ))
  if(nrow(utr5Hits)>0){
    utr5Hits <- cbind(snpInfo[utr5Hits$QUERYID,], utr5Hits[,.(tx_id=as.integer(TXID),type="utr5")])
    utr5Hits <- merge(utr5Hits, txinfo[,.(tx_id,anno=tx_name)],by="tx_id")[,-c("tx_id")]
    # collapse annotation:
    utr5Hits <- utr5Hits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  }else{
    utr5Hits <- cpgHits[0,]
  }
  message("      5'utr hits: ",nrow(utr5Hits))

  # annotate snp with 3 utr:
  utr3Hits <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges,
                                                                                                 txdb,
                                                                                                 VariantAnnotation::ThreeUTRVariants())) ))
  if(nrow(utr3Hits)>0){
    utr3Hits <- cbind(snpInfo[utr3Hits$QUERYID,], utr3Hits[,.(tx_id=as.integer(TXID),type="utr3")])
    utr3Hits <- merge(utr3Hits, txinfo[,.(tx_id,anno=tx_name)],by="tx_id")[,-c("tx_id")]
    # collapse annotation:
    utr3Hits <- utr3Hits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  }else{
    utr3Hits <- cpgHits[0,]
  }
  message("      3'utr hits: ",nrow(utr3Hits))

  # annotate snp with splice sites:
  splicingHits <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges,
                                                                                                     txdb,
                                                                                                     VariantAnnotation::SpliceSiteVariants())) ))
  if(nrow(splicingHits)>0){
    splicingHits <- cbind(snpInfo[splicingHits$QUERYID,], splicingHits[,.(tx_id=as.integer(TXID),type="spliceSite")])
    splicingHits <- merge(splicingHits, txinfo[,.(tx_id,anno=tx_name)],by="tx_id")[,-c("tx_id")]
    # collapse annotation:
    splicingHits <- splicingHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  }else{
    splicingHits <- cpgHits[0,]
  }
  message("      splice site hits: ",nrow(splicingHits))

  intergenicHits <- suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges,
                                                                                      txdb,
                                                                                      VariantAnnotation::IntergenicVariants(), ignore.strand=TRUE)) )
  if(nrow(intergenicHits)>0){
    fetchLastElement <- function(x){ a=unlist(as.list(x)); ifelse(length(a)==0,"", a[length(a)]) }
    fetchFirstElement <- function(x){ a=unlist(as.list(x));  ifelse(length(a)==0,"", a[1]) }
    intergenicHits[,c("anno", "type"):=.(paste0(unlist(lapply(intergenicHits$PRECEDEID, fetchLastElement)),"-",unlist(lapply(intergenicHits$FOLLOWID, fetchFirstElement))), "intergenic")]
    intergenicHits <- cbind(snpInfo[intergenicHits$QUERYID,], intergenicHits[,.(anno,type)])
    # collapse annotation:
    intergenicHits <- intergenicHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
    rm(fetchLastElement, fetchFirstElement)
  }else{
    intergenicHits <- cpgHits[0,]
  }
  message("      Intergenic hits: ",nrow(promoterHits))
  rm(txinfo)

  # combine all hits:
  snpHits <- do.call(rbind,list(cpgHits, enhancerHits, promoterHits, exonHits, cdsHits, intronHits, utr3Hits, utr5Hits, tfHits, splicingHits, intergenicHits))
  # For the rest of the failed to be annotated, all in intergenic:
  failedSnps <- fsetdiff(snpInfo, snpHits[,.(chrom, pos, pValue)])
  failedSnps$anno <- "-"
  failedSnps$type <- "intergenic"
  snpHits <- rbind(snpHits, failedSnps)

  # prop:
  hitsProp <- as.data.table(round(prop.table(table(snpHits$type))*100, 3))
  names(hitsProp) <- c("type", "proportion")
  hitsProp <- hitsProp[order(-proportion)]
  message("==> Proportion: \n", paste("     ",paste(hitsProp$type, ":\t",hitsProp$proportion, "%",sep=""), collapse = "; \n"))

  return(snpHits)
}



#' @title Compare P-values reported in the association result file to P-values calculated from Z statistics derived from the reported beta and standard error.
#'
#' @param summaryDT A data.frame with three cols: pval,  beta, se.
#'
#' @return a list containing a data.frame of estimated pvalues and A ggplot2 object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "https://raw.githubusercontent.com/dingruofan/exampleData/master/eqtl/MMP7_qtlDF.txt"
#' qtl <- fread(url1, sep="\t")
#' xQTLanno_PZPlot(qtl[,.(pValue, beta, se)])
#' }
xQTLanno_PZPlot <- function(summaryDT){
  se <- pval <- pZtest <- logPZ <- logP <- NULL
  . <-NULL
  summaryDT <- na.omit(summaryDT)
  summaryDT <- summaryDT[,1:3]
  names(summaryDT) <- c("pval", "beta", "se")
  summaryDT[,c("pZtest", "logP") := .(pnorm(-abs(beta/se)), log(pval, 10)*(-1))]
  summaryDT[,"logPZ" := log(pZtest, 10)*(-1)]
  summaryDT <- na.omit(summaryDT)

  p <- ggplot(summaryDT)+
    geom_point(aes(x=logPZ, y=logP))+
    geom_abline(intercept = 0)+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))+
    ylab(expression(-log["10"]("Pvalue-raw")))+
    xlab(expression(-log["10"]("Pvalue-estimated")))+
    theme_classic()+
    theme(
      axis.text = element_text(rel(1.3)),
      axis.title = element_text(rel(1.4))
    )
  print(p)
  return(list(summaryDT=summaryDT, p=p))
}
