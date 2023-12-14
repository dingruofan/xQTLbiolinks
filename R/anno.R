#' @title Calculate genomic control inflation factor for a QTL/GWAS summary statistics dataset
#'
#' @param summaryDT A data.frame containing one or two columns: p-value (required) and group (optional)
#' @import data.table
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/eqtl/MMP7_qtlDF.txt"
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


#' @title Annotate GWAS / QTL signals using bed4 format table object that provided by user
#'
#' @param snpInfo A data.table/data.frame with three columns: chromosome, position and p-value.
#' @param genomeVersion "hg38" (default) or "hg19". Note: hg19 will be converted to hg38 automatically.
#' @param enrichElement A data.table of data.frame object including 4 columns (consistent with bed4 format): chrom, start, end, name.
#' @param distLimit Defaults: 0 (variants overlap with elements).
#' @importFrom GenomicRanges distanceToNearest mcols
#' @import data.table
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwas/gwasSub.txt.gz"
#' url2 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/enhancer.txt"
#' snpInfo <- data.table::fread(url1, sep="\t")
#' enhancerDT <- data.table::fread(url2, sep="\t")
#' variants_hit_enhancer <- xQTLanno_chippeak(snpInfo, enrichElement=enhancerDT)
#' }
xQTLanno_chippeak <- function(snpInfo="",  genomeVersion="hg38", enrichElement=NULL, distLimit=1){
  chrom <- chr <- chromStart <- chromEnd <- elementType <- pValue <- tx_name <- distance <- V1 <- V2 <- V3 <- NULL
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

  #
  message("==> Start annotating...")
  data.table::setDT(enrichElement)
  names(enrichElement) <- c("chr", "chromStart", "chromEnd", "elementType")
  enrichElement <- enrichElement[,.(chr= as.character(chr), chromStart = as.numeric(chromStart),chromEnd = as.numeric(chromEnd), elementType=as.character(elementType))]
  enrichElementRange <- GenomicRanges::GRanges(enrichElement$chr,
                                               IRanges::IRanges(enrichElement$chromStart, enrichElement$chromEnd),
                                               strand= "*",
                                               enrichElement[,"elementType"]
  )
  nearestDist <- data.table::as.data.table(GenomicRanges::distanceToNearest(snpRanges, enrichElementRange))
  if(nrow(nearestDist)>0){
    nearestDist <- do.call(cbind, list(snpInfo[nearestDist$queryHits,],  data.table::as.data.table(GenomicRanges::mcols(enrichElementRange))[nearestDist$subjectHits,.(anno=elementType)],nearestDist[,.(dist=distance)] ))
    nearestDist <- nearestDist[dist<=distLimit,][order(dist)]
    return(nearestDist)
  }
  else{
    return(NULL)
  }
}


#' @title Genomic annotation of significant signals from GWAS / QTL
#'
#' @param snpInfo A data.table/data.frame with three columns: chromosome, position, and P-value.
#' @param genomeVersion "hg38" (default) or "hg19". Note: hg19 will be converted to hg38 automatically.
#' @param p_cutoff Cutoff of p-values of significant variants that will be annotated
#' @import data.table
#' @import stringr
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom GenomicFeatures intronicParts exonicParts transcripts
#' @importFrom IRanges findOverlaps IRanges
#'
#' @return A data.table object of variants' genomics distribution
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwas/gwasSub.txt.gz"
#' snpInfo <- data.table::fread(url1, sep="\t")
#' snpHits <- xQTLanno_genomic(snpInfo)
#' }
xQTLanno_genomic <- function(snpInfo="", p_cutoff =5e-8, genomeVersion="hg38"){

  width <- N <- enrichment <- p.value <- NULL
  chrom <- pValue <- V1 <- V2 <- V3 <- tx_name <- anno <- TXID <- tx_id <- type <- pos <- proportion <- NULL
  nrow_id <- rep_id<- NULL

  datatype_DT <- data.table(type=c("spliceSite",  "utr5","utr3", "exon", "promoter", "intron", "downstream", "intergenic"),
                            bg_ratio = c(0.017,    0.3,   0.98,   2.13,   1.09,       50.27,    1.1,          44.11)/100)

  .<- NULL

  if( !requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene") ){ stop("please install package \"TxDb.Hsapiens.UCSC.hg38.knownGene\" from BioConductor.") }

  message("==> Start checking variants...")
  # check snpinfo:
  snpInfo <- data.table::as.data.table(snpInfo[,1:3])
  snpInfo <- na.omit(snpInfo)
  names(snpInfo) <- c("chrom", "pos", "pValue")
  snpInfo$chrom <- as.character(snpInfo$chrom)
  if(!stringr::str_detect(snpInfo[1]$chrom, "chr")){snpInfo$chrom <- paste0("chr", snpInfo$chrom)}
  snpInfo$pos <- as.integer(snpInfo$pos)
  snpInfo$pValue <- as.numeric(snpInfo$pValue)
  snpInfo <- snpInfo[chrom %in% paste0("chr",c(1:22,"X","Y","M"))]
  snpInfo <- snpInfo[pValue>0,]
  if(nrow(snpInfo)<1){
    stop("Number of variants < 1...")
  }
  snpInfo_non_sig <- snpInfo[pValue >= p_cutoff,]
  snpInfo <- snpInfo[pValue<p_cutoff,]
  message("    Number of SNPs with p-value <", p_cutoff,": ",nrow(snpInfo))
  # create snpinfo Grange object:
  snpRanges <- GenomicRanges::GRanges(snpInfo$chrom,
                                      IRanges::IRanges(snpInfo$pos, snpInfo$pos),
                                      strand= rep("*", nrow(snpInfo)),
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

  # all transcripts:
  txinfo <- as.data.table(GenomicFeatures::transcripts(txdb))
  txinfo <- txinfo[which(seqnames %in% paste0("chr", c(1:22,"X","Y", "M"))),]

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
  promoterDT <- as.data.table(promoterParts)
  promoterDT$TSS <- ifelse(promoterDT$strand=="+", promoterDT$start+2000, promoterDT$end-2000)
  promoterDT$upstream <- ifelse(promoterDT$strand=="+", promoterDT$TSS-1000, promoterDT$TSS+1000)
  promoterDT$start <- ifelse(promoterDT$strand=="+", promoterDT$upstream, promoterDT$TSS)
  promoterDT$end <- ifelse(promoterDT$strand=="+", promoterDT$TSS, promoterDT$upstream)
  promoterParts <- GenomicRanges::GRanges(as.character(promoterDT$seqnames),
                                          IRanges::IRanges(promoterDT$start, promoterDT$end),
                                          strand= promoterDT$strand,
                                          promoterDT[,.(tx_name)]
  )

  # extract downstream:
  txinfo$downstream <- ifelse(txinfo$strand=="+", txinfo$end+1000, txinfo$start-1000)
  txinfo$downstream_start <- ifelse(txinfo$strand=="+", txinfo$end, txinfo$downstream)
  txinfo$downstream_end <- ifelse(txinfo$strand=="+", txinfo$downstream, txinfo$start)

  downstreamParts <- GenomicRanges::GRanges(as.character(txinfo$seqnames),
                                          IRanges::IRanges(txinfo$downstream_start, txinfo$downstream_end),
                                          strand= txinfo$strand,
                                          txinfo[,.(tx_name)]
  )


  #####################
  # construct non_sig region:
  snpRanges_non_sig_all <- GenomicRanges::GRanges(snpInfo_non_sig$chrom,
                                                  IRanges::IRanges(snpInfo_non_sig$pos, snpInfo_non_sig$pos),
                                                  strand= rep("*", nrow(snpInfo_non_sig)),
                                                  snpInfo_non_sig[,.(pValue)]
  )

  ############ start find overlap with annotations:
  message("==> Start annotating variants...")

  # default:
  cpgHits <- data.table(chrom=character(0),pos=numeric(0),pValue=numeric(0),type=character(0),anno=character(0))

  # annotate snp with introns:
  intronHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                                 subject = intronParts,
                                                                 maxgap = 1,
                                                                 ignore.strand=TRUE ))
  convertCharacter <- function(x){ a=as.list(x);paste0(unlist(a), collapse = ",") }
  if(nrow(intronHits)>0){
    intronHits <- cbind(snpInfo[intronHits$queryHits,],intronAnno[intronHits$subjectHits,.(anno=tx_name, type="intron")])
    intronHits$anno <- unlist(lapply(intronHits$anno, convertCharacter))
    # collapse annotation:
    intronHits <- intronHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  }else{
    intronHits <- cpgHits[0,]
  }
  intron_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig_all,
                                                                          subject = intronParts,
                                                                          maxgap = 1,
                                                                          ignore.strand=TRUE ))
  intron_hits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[intron_hits_non_sig_all$queryHits]))
  intron_hits_non_sig_all$type<- "intron"
  # message("      intron hits: ",nrow(intronHits),"; none-hits: ",nrow(intron_hits_non_sig_all))
  # rm(intronParts, intronAnno)

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
  exon_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig_all,
                                                                        subject = exonParts,
                                                                        maxgap = 1,
                                                                        ignore.strand=TRUE ))
  exon_hits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[exon_hits_non_sig_all$queryHits]))
  exon_hits_non_sig_all$type <- "exon"
  # message("      exon hits: ",nrow(exonHits), "; none-hits: ", nrow(exon_hits_non_sig_all))
  # rm(exonParts, exonAnno)

  # annotate snp with promoters:
  promoterHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                                   subject = promoterParts,
                                                                   maxgap = 1,
                                                                   ignore.strand=TRUE ))
  promoterHits <- cbind(snpInfo[promoterHits$queryHits,], as.data.table(GenomicRanges::mcols(promoterParts))[promoterHits$subjectHits,.(anno=tx_name, type="promoter")])
  # collapse annotation:
  promoterHits <- promoterHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  promoter_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig_all,
                                                                            subject = promoterParts,
                                                                            maxgap = 1,
                                                                            ignore.strand=TRUE ))
  promoter_hits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[promoter_hits_non_sig_all$queryHits]))
  promoter_hits_non_sig_all$type <- "promoter"


  # message("      promoter hits: ",nrow(promoterHits), "; none-hits: ", nrow(promoter_hits_non_sig_all))
  # rm(promoterParts)

  # annotate snp with downstream:
  downstreamHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                                   subject = downstreamParts,
                                                                   maxgap = 1,
                                                                   ignore.strand=TRUE ))
  downstreamHits <- cbind(snpInfo[downstreamHits$queryHits,], as.data.table(GenomicRanges::mcols(downstreamParts))[downstreamHits$subjectHits,.(anno=tx_name, type="downstream")])
  # collapse annotation:
  downstreamHits <- downstreamHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  downstream_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig_all,
                                                                                subject = downstreamParts,
                                                                                maxgap = 1,
                                                                                ignore.strand=TRUE ))
  downstream_hits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[downstream_hits_non_sig_all$queryHits]))
  downstream_hits_non_sig_all$type <- "downstream"

  # message("      downstream hits: ",nrow(downstreamHits), "; none-hits: ", nrow(downstream_hits_non_sig_all))
  # rm(promoterParts)


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
  utr5_hits_non_sig_all <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig_all,
                                                                                                 txdb,
                                                                                                 VariantAnnotation::FiveUTRVariants())) ))
  utr5_hits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[utr5_hits_non_sig_all$QUERYID]))
  utr5_hits_non_sig_all$type <- "utr5"
  # message("      5'utr hits: ",nrow(utr5Hits),"; none hits: ", nrow(utr5_hits_non_sig_all))

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
  utr3_hits_non_sig_all <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig_all,
                                                                                                          txdb,
                                                                                                          VariantAnnotation::ThreeUTRVariants())) ))
  utr3_hits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[utr3_hits_non_sig_all$QUERYID]))
  utr3_hits_non_sig_all$type <- "utr3"
  # message("      3'utr hits: ",nrow(utr3Hits),"; none-hits: ", nrow(utr3_hits_non_sig_all))

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
  splicing_hits_non_sig_all <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig_all,
                                                                                                              txdb,
                                                                                                              VariantAnnotation::SpliceSiteVariants())) ))
  splicing_hits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[splicing_hits_non_sig_all$QUERYID]))
  splicing_hits_non_sig_all$type <- "spliceSite"
  # message("      splice site hits: ",nrow(splicingHits),"; none-hits: ", nrow(splicing_hits_non_sig_all))

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
  intergenicHits_non_sig_all <- suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig_all,
                                                                                      txdb,
                                                                                      VariantAnnotation::IntergenicVariants(), ignore.strand=TRUE)) )
  intergenicHits_non_sig_all <- unique(as.data.table(snpRanges_non_sig_all[intergenicHits_non_sig_all$QUERYID]))
  intergenicHits_non_sig_all$type <- "intergenic"
  # message("      Intergenic hits: ",nrow(promoterHits), "; none-hits: ", nrow(intergenicHits_non_sig_all))
  # rm(txinfo)

  # combine all hits:
  snpHits <- do.call(rbind,list(promoterHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                downstreamHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                exonHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                # cdsHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                intronHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                utr3Hits[,c("chrom", "pos", "pValue", "anno", "type")],
                                utr5Hits[,c("chrom", "pos", "pValue", "anno", "type")],
                                splicingHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                intergenicHits[,c("chrom", "pos", "pValue", "anno", "type")]))
  # For the rest of the failed to be annotated, all in intergenic:
  failedSnps <- fsetdiff(snpInfo, snpHits[,.(chrom, pos, pValue)])
  failedSnps$anno <- "-"
  failedSnps$type <- "intergenic"
  snpHits <- rbind(snpHits, failedSnps)

  snpHits <- merge(datatype_DT, snpHits, by="type", sort=FALSE)
  snpHits <- snpHits[,.SD[1,], by=c("chrom", "pos", "pValue")]
  # prop:
  hitsProp <- round(prop.table(table(snpHits$type))*100, 3)
  hitsProp <- data.table(type=names(hitsProp), proportion=as.numeric(hitsProp))
  hitsProp <- hitsProp[order(-proportion)]
  message("==> Proportion: \n", paste("     ",paste(hitsProp$type, ":\t",hitsProp$proportion, "%",sep=""), collapse = "; \n"))
  # enrichment:
  hitsProp <- merge(datatype_DT, hitsProp, by="type")
  hitsProp$enrichment <- hitsProp$proportion/100 / hitsProp$bg_ratio


  # non-sig variants located in element(no random):
  non_sig_all <- rbind(intron_hits_non_sig_all, downstream_hits_non_sig_all, exon_hits_non_sig_all, promoter_hits_non_sig_all, utr5_hits_non_sig_all, utr3_hits_non_sig_all, splicing_hits_non_sig_all, intergenicHits_non_sig_all)
  dup_non_sig <- unique(non_sig_all[duplicated(non_sig_all[,.(seqnames, start, pValue)]),.(seqnames, start, pValue)])
  dup_non_sig <- merge(non_sig_all, dup_non_sig, by=c("seqnames", "start", "pValue"), sort=FALSE)[,.(seqnames,start,end,width,strand,pValue,type)]
  noDup_non_sig <- fsetdiff(non_sig_all, dup_non_sig)
  dup_non_sig <- merge(datatype_DT, dup_non_sig, by="type", sort=FALSE)
  dup_non_sig <- dup_non_sig[,.SD[1,], by=c("seqnames", "start", "pValue")]
  non_sig_all <- rbind(dup_non_sig[,-c("bg_ratio")], noDup_non_sig)

  non_sig_count <- as.data.table(table(non_sig_all$type))
  names(non_sig_count) <- c("type", "num_non_sig_hits")

  # non_sig_count <- data.table(
  #   type = c("intron", "downstream", "exon", "promoter",  "utr5", "utr3", "spliceSite", "intergenic"),
  #   num_non_sig_hits = c(nrow(intron_hits_non_sig_all), nrow(downstream_hits_non_sig_all), nrow(exon_hits_non_sig_all), nrow(promoter_hits_non_sig_all), nrow(utr5_hits_non_sig_all), nrow(utr3_hits_non_sig_all), nrow(splicing_hits_non_sig_all), nrow(intergenicHits_non_sig_all))
  #   )

  snpHits_count  <- merge(as.data.table(table(snpHits$type))[,.(type=V1, num_sig_hits=N)], non_sig_count, by = "type")
  snpHits_count$num_sig_non_hits <- nrow(snpInfo) - snpHits_count$num_sig_hits
  snpHits_count$num_non_sig_non_hits <- nrow(snpInfo_non_sig) - snpHits_count$num_non_sig_hits
  snpHits_count$p.value = apply(snpHits_count, 1, function(x){ x<-as.numeric(x[2:5]); x_enrich = fisher.test(matrix(c(x[1],x[2],x[3],x[4]), nrow=2, byrow = TRUE),alternative = "two.sided"); return(x_enrich$p.value)})

  snpEnrich <- merge(hitsProp[,.(type, enrichment)], snpHits_count[,.(type, p.value)], by="type", sort=FALSE)[order(enrichment)]

  return(list(snpHits = snpHits, snpEnrich = snpEnrich, snpHits_count=snpHits_count[,-c("p.value")]))
}

# snpHits <- xQTLanno_genomic(gwasDT[,.(chr, position,`p-value`)], p_cutoff = 5e-8, genomeVersion = "hg19")
# annoPlotBar <- xQTLvisual_anno(snpHits2,plotType="bar")






