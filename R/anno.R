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
#' url1 <- "https://master.dl.sourceforge.net/project/exampledata/eqtl/MMP7_qtlDF.txt"
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
#' @param snpInfo A data.table/data.frame with three columns: chromosome, position and p-value.
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
#' url1 <- "https://master.dl.sourceforge.net/project/exampledata/gwas/gwasSub.txt.gz"
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
    # tfAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "tf_filtered_sorted_merge_sub.bed.gz"), header=FALSE)
    temp1 <- tempfile(fileext=".gz")
    download.file("https://master.dl.sourceforge.net/project/exampledata/exdata/tf_filtered_sorted_merge_sub.bed.gz", temp1)
    close(file(temp1))
    tfAnno <- fread(temp1, header=FALSE)
    rm(temp1)
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
#' url1 <- "https://master.dl.sourceforge.net/project/exampledata/gwas/gwasSub.txt.gz"
#' snpInfo <- data.table::fread(url1, sep="\t")
#' snpHits <- xQTLanno_genomic(snpInfo)
#' }
xQTLanno_genomic <- function(snpInfo="", p_cutoff =5e-8, genomeVersion="hg38"){
  random_repeat_times <- 500
  chrom <- pValue <- V1 <- V2 <- V3 <- tx_name <- anno <- TXID <- tx_id <- type <- pos <- proportion <- NULL
  .<- NULL
  message("==> Start checking variants...")
  # check snpinfo:
  snpInfo <- data.table::as.data.table(snpInfo[,1:3])
  snpInfo <- na.omit(snpInfo)
  names(snpInfo) <- c("chrom", "pos", "pValue")
  snpInfo$chrom <- as.character(snpInfo$chrom)
  if(!str_detect(snpInfo[1]$chrom, "chr")){snpInfo$chrom <- paste0("chr", snpInfo$chrom)}
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
  message("      cpg...")
  # cpgAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "cpgLsland.bed.gz"), header=FALSE)
  # temp1 <- tempfile(fileext=".gz")
  # download.file("https://master.dl.sourceforge.net/project/exampledata/cpgLsland.bed.gz?viasf=1", temp1)
  # close(file(temp1))
  # cpgAnno <- fread(temp1, header=FALSE)
  # rm(temp1)
  cpgAnno <- fread("https://master.dl.sourceforge.net/project/exampledata/cpgLsland.bed.gz", header=FALSE)


  cpgRanges <- GenomicRanges::GRanges(cpgAnno$V1,
                                      IRanges::IRanges(cpgAnno$V2, cpgAnno$V3), strand= "*")

  # import enhancer dataset:
  message("      enhancer...")
  # enhancerAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "geneHancer.bed.gz"), header=FALSE)
  # temp1 <- tempfile(fileext=".gz")
  # # download.file("https://master.dl.sourceforge.net/project/exampledata/exdata/geneHancer.bed.gz", temp1)
  # download.file("https://master.dl.sourceforge.net/project/exampledata/geneHancer.bed.gz?viasf=1", temp1)
  # close(file(temp1))
  # enhancerAnno <- fread(temp1, header=FALSE)
  # rm(temp1)
  enhancerAnno <- fread("https://master.dl.sourceforge.net/project/exampledata/geneHancer.bed.gz", header = FALSE)
  enhancerRanges <- GenomicRanges::GRanges(enhancerAnno$V1,
                                           IRanges::IRanges(enhancerAnno$V2, enhancerAnno$V3), strand= "*")

  # import tf cluster dataset:
  message("      tf cluster...")
  # tfAnno <- fread(system.file(package = "xQTLbiolinks", "extdata", "tf_filtered_sorted_merge_sub.bed.gz"), header=FALSE)
  # temp1 <- tempfile(fileext=".gz")
  # # download.file("https://github.com/dingruofan/exampleData/raw/master/exdata/tf_filtered_sorted_merge_sub.bed.gz", temp1)
  # download.file("https://master.dl.sourceforge.net/project/exampledata/tf_filtered_sorted_merge_sub.bed.gz?viasf=1", temp1)
  # close(file(temp1))
  # tfAnno <- fread(temp1, header=FALSE)
  tfAnno <- fread("https://master.dl.sourceforge.net/project/exampledata/tf_filtered_sorted_merge_sub.bed.gz", header=FALSE)
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
  txinfo <- as.data.table(GenomicFeatures::transcripts(txdb))
  txinfo <- txinfo[which(seqnames %in% paste0("chr", c(1:22,"X","Y", "M"))),]

  #####################
  # construct random signals:
  non_sig_id_random <- rbindlist(lapply(1:random_repeat_times, function(x){data.table(rep_id = x, nrow_id = sample(1:nrow(snpInfo_non_sig), nrow(snpInfo)))}))
  non_sig_id_random_uniq <- cbind(snpInfo_non_sig[unique(non_sig_id_random$nrow_id),],unique(non_sig_id_random[,.(nrow_id)]))
  snpRanges_non_sig <- GenomicRanges::GRanges(non_sig_id_random_uniq$chrom,
                                              IRanges::IRanges(non_sig_id_random_uniq$pos, non_sig_id_random_uniq$pos),
                                              strand= rep("*", nrow(non_sig_id_random_uniq)),
                                              non_sig_id_random_uniq[,.(pValue, nrow_id)]
  )

  ############ start find overlap with annotations:
  message("==> Start annotating variants...")
  cpgHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                              subject = cpgRanges,
                                                              maxgap = 1,
                                                              ignore.strand=TRUE ))
  cpgHits <- cbind(snpInfo[cpgHits$queryHits,],cpgAnno[cpgHits$subjectHits,.(anno=paste0(V1,":",V2,"-",V3), type="cpg")])
  cpg_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig,
                                                                       subject = cpgRanges,
                                                                       maxgap = 1,
                                                                       ignore.strand=TRUE ))
  cpg_hits_non_sig_all <- as.data.table(snpRanges_non_sig[cpg_hits_non_sig_all$queryHits])
  cpg_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    cpg_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], cpg_hits_non_sig_all, by="nrow_id"))
    non_cpg_hits <- nrow(snpInfo) - nrow(cpgHits)
    non_cpg_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - cpg_hits_non_sig
    cpg_confusion <- matrix(c(nrow(cpgHits), non_cpg_hits, cpg_hits_non_sig, non_cpg_hits_non_sig), nrow=2)
    cpg_enrich <- fisher.test(cpg_confusion, alternative = "greater")
    cpg_enrich <- data.table( type="cpg",p.value = cpg_enrich$p.value,  OR=cpg_enrich$estimate)
    return(cpg_enrich)
  }))
  message("      cpg hits: ",nrow(cpgHits))
  rm(cpgAnno, cpgRanges)

  # annotate snp with enhancers:
  enhancerHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                                   subject = enhancerRanges,
                                                                   maxgap = 1,
                                                                   ignore.strand=TRUE ))
  enhancerHits <- cbind(snpInfo[enhancerHits$queryHits,],enhancerAnno[enhancerHits$subjectHits,.(anno=paste0(V1,":",V2,"-",V3), type="enhancer")])
  enhancer_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig,
                                                                           subject = enhancerRanges,
                                                                           maxgap = 1,
                                                                           ignore.strand=TRUE ))
  enhancer_hits_non_sig_all <- as.data.table(snpRanges_non_sig[enhancer_hits_non_sig_all$queryHits])
  enhancer_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    enhancer_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], enhancer_hits_non_sig_all, by="nrow_id"))
    non_enhancer_hits <- nrow(snpInfo) - nrow(enhancerHits)
    non_enhancer_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - enhancer_hits_non_sig
    enhancer_confusion <- matrix(c(nrow(enhancerHits), non_enhancer_hits, enhancer_hits_non_sig, non_enhancer_hits_non_sig), nrow=2)
    enhancer_enrich <- fisher.test(enhancer_confusion, alternative = "greater")
    enhancer_enrich <- data.table( type="enhancer",p.value = ifelse(enhancer_enrich$p.value==0,2.2e-16, enhancer_enrich$p.value),  OR=enhancer_enrich$estimate)
    return(enhancer_enrich)
  }))
  message("      enhancer hits: ",nrow(enhancerHits))
  rm(enhancerAnno, enhancerRanges)

  # annotate snp with tf clusters:
  tfHits <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges,
                                                             subject = tfRanges,
                                                             maxgap = 1,
                                                             ignore.strand=TRUE ))
  tfHits <- cbind(snpInfo[tfHits$queryHits,],tfAnno[tfHits$subjectHits,.(anno=paste0(V1,":",V2,"-",V3), type="tfCluster")])
  tf_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig,
                                                                          subject = tfRanges,
                                                                          maxgap = 1,
                                                                          ignore.strand=TRUE ))
  tf_hits_non_sig_all <- as.data.table(snpRanges_non_sig[tf_hits_non_sig_all$queryHits])
  tf_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    tf_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], tf_hits_non_sig_all, by="nrow_id"))
    non_tf_hits <- nrow(snpInfo) - nrow(tfHits)
    non_tf_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - tf_hits_non_sig
    tf_confusion <- matrix(c(nrow(tfHits), non_tf_hits, tf_hits_non_sig, non_tf_hits_non_sig), nrow=2)
    tf_enrich <- fisher.test(tf_confusion, alternative = "greater")
    tf_enrich <- data.table( type="tfCluster",p.value = ifelse(tf_enrich$p.value==0,2.2e-16, tf_enrich$p.value),  OR=tf_enrich$estimate)
    return(tf_enrich)
  }))

  message("      tf hits: ",nrow(tfHits))
  rm(tfAnno, tfRanges)

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
  intron_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig,
                                                                          subject = intronParts,
                                                                          maxgap = 1,
                                                                          ignore.strand=TRUE ))
  intron_hits_non_sig_all <- as.data.table(snpRanges_non_sig[intron_hits_non_sig_all$queryHits])
  intron_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    intron_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], intron_hits_non_sig_all, by="nrow_id", sort=FALSE))
    non_intron_hits <- nrow(snpInfo) - nrow(intronHits)
    non_intron_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - intron_hits_non_sig
    intron_confusion <- matrix(c(nrow(intronHits), non_intron_hits, intron_hits_non_sig, non_intron_hits_non_sig), nrow=2)
    intron_enrich <- fisher.test(intron_confusion, alternative = "greater")
    intron_enrich <- data.table( type="intron",p.value = ifelse(intron_enrich$p.value==0,2.2e-16, intron_enrich$p.value),  OR=intron_enrich$estimate)
    return(intron_enrich)
  }))
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
  exon_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig,
                                                                        subject = exonParts,
                                                                        maxgap = 1,
                                                                        ignore.strand=TRUE ))
  exon_hits_non_sig_all <- as.data.table(snpRanges_non_sig[exon_hits_non_sig_all$queryHits])
  exon_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    exon_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], exon_hits_non_sig_all, by="nrow_id"))
    non_exon_hits <- nrow(snpInfo) - nrow(exonHits)
    non_exon_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - exon_hits_non_sig
    exon_confusion <- matrix(c(nrow(exonHits), non_exon_hits, exon_hits_non_sig, non_exon_hits_non_sig), nrow=2)
    exon_enrich <- fisher.test(exon_confusion, alternative = "greater")
    exon_enrich <- data.table( type="exon",p.value = ifelse(exon_enrich$p.value==0,2.2e-16, exon_enrich$p.value),  OR=exon_enrich$estimate)
    return(exon_enrich)
  }))

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
  promoter_hits_non_sig_all <- as.data.table(SummarizedExperiment::findOverlaps(query = snpRanges_non_sig,
                                                                            subject = promoterParts,
                                                                            maxgap = 1,
                                                                            ignore.strand=TRUE ))
  promoter_hits_non_sig_all <- as.data.table(snpRanges_non_sig[promoter_hits_non_sig_all$queryHits])
  promoter_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    promoter_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], promoter_hits_non_sig_all, by="nrow_id"))
    non_promoter_hits <- nrow(snpInfo) - nrow(promoterHits)
    non_promoter_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - promoter_hits_non_sig
    promoter_confusion <- matrix(c(nrow(promoterHits), non_promoter_hits, promoter_hits_non_sig, non_promoter_hits_non_sig), nrow=2)
    promoter_enrich <- fisher.test(promoter_confusion, alternative = "greater")
    promoter_enrich <- data.table( type="promoter",p.value = ifelse(promoter_enrich$p.value==0,2.2e-16, promoter_enrich$p.value),  OR=promoter_enrich$estimate)
    return(promoter_enrich)
  }))

  message("      promoter hits: ",nrow(promoterHits))
  rm(promoterParts)

  # annotate snp with cds:
  # cdsHits <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges,
  #                                                                                               txdb,
  #                                                                                               VariantAnnotation::CodingVariants())) ))
  # if(nrow(cdsHits)>0){
  #   cdsHits <- cbind(snpInfo[cdsHits$QUERYID,], cdsHits[,.(tx_id=as.integer(TXID),type="nonsynonymous")])
  #   cdsHits <- merge(cdsHits, txinfo[,.(tx_id,anno=tx_name)],by="tx_id")[,-c("tx_id")]
  #   # collapse annotation:
  #   cdsHits <- cdsHits[,.(anno=paste(anno,collapse = ",")),by=c("chrom", "pos", "pValue", "type")]
  # }else{
  #   cdsHits <- cpgHits[0,]
  # }
  # message("      nonsynonymous hits: ",nrow(cdsHits))

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
  utr5_hits_non_sig_all <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig,
                                                                                                 txdb,
                                                                                                 VariantAnnotation::FiveUTRVariants())) ))
  utr5_hits_non_sig_all <- as.data.table(snpRanges_non_sig[utr5_hits_non_sig_all$QUERYID])
  utr5_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    utr5_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], utr5_hits_non_sig_all, by="nrow_id"))
    non_utr5_hits <- nrow(snpInfo) - nrow(utr5Hits)
    non_utr5_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - utr5_hits_non_sig
    utr5_confusion <- matrix(c(nrow(utr5Hits), non_utr5_hits, utr5_hits_non_sig, non_utr5_hits_non_sig), nrow=2)
    utr5_enrich <- fisher.test(utr5_confusion, alternative = "greater")
    utr5_enrich <- data.table( type="utr5",p.value = ifelse(utr5_enrich$p.value==0,2.2e-16, utr5_enrich$p.value),  OR=utr5_enrich$estimate)
    return(utr5_enrich)
  }))
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
  utr3_hits_non_sig_all <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig,
                                                                                                          txdb,
                                                                                                          VariantAnnotation::ThreeUTRVariants())) ))
  utr3_hits_non_sig_all <- as.data.table(snpRanges_non_sig[utr3_hits_non_sig_all$QUERYID])

  utr3_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    utr3_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], utr3_hits_non_sig_all, by="nrow_id"))
    non_utr3_hits <- nrow(snpInfo) - nrow(utr3Hits)
    non_utr3_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - utr3_hits_non_sig
    utr3_confusion <- matrix(c(nrow(utr3Hits), non_utr3_hits, utr3_hits_non_sig, non_utr3_hits_non_sig), nrow=2)
    utr3_enrich <- fisher.test(utr3_confusion, alternative = "greater")
    utr3_enrich <- data.table( type="utr3",p.value = ifelse(utr3_enrich$p.value==0,2.2e-16, utr3_enrich$p.value),  OR=utr3_enrich$estimate)
    return(utr3_enrich)
  }))

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
  splicing_hits_non_sig_all <- suppressMessages(suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig,
                                                                                                              txdb,
                                                                                                              VariantAnnotation::SpliceSiteVariants())) ))
  splicing_hits_non_sig_all <- as.data.table(snpRanges_non_sig[splicing_hits_non_sig_all$QUERYID])

  splicing_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    splicing_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], splicing_hits_non_sig_all, by="nrow_id"))
    non_splicing_hits <- nrow(snpInfo) - nrow(splicingHits)
    non_splicing_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - splicing_hits_non_sig
    splicing_confusion <- matrix(c(nrow(splicingHits), non_splicing_hits, splicing_hits_non_sig, non_splicing_hits_non_sig), nrow=2)
    splicing_enrich <- fisher.test(splicing_confusion, alternative = "greater")
    splicing_enrich <- data.table( type="splicing",p.value = ifelse(splicing_enrich$p.value==0,2.2e-16, splicing_enrich$p.value),  OR=splicing_enrich$estimate)
    return(splicing_enrich)
  }))
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
  intergenicHits_non_sig_all <- suppressWarnings( as.data.table(VariantAnnotation::locateVariants(snpRanges_non_sig,
                                                                                      txdb,
                                                                                      VariantAnnotation::IntergenicVariants(), ignore.strand=TRUE)) )
  intergenic_hits_non_sig_all <- as.data.table(snpRanges_non_sig[intergenicHits_non_sig_all$QUERYID])
  intergenic_enrich <- rbindlist(lapply(1:random_repeat_times, function(x){
    # cat("| ",x)
    intergenic_hits_non_sig <- nrow(merge(non_sig_id_random[rep_id == x,], intergenic_hits_non_sig_all, by="nrow_id"))
    non_intergenic_hits <- nrow(snpInfo) - nrow(intergenicHits)
    non_intergenic_hits_non_sig <- nrow(non_sig_id_random[rep_id == x,]) - intergenic_hits_non_sig
    intergenic_confusion <- matrix(c(nrow(intergenicHits), non_intergenic_hits, intergenic_hits_non_sig, non_intergenic_hits_non_sig), nrow=2)
    intergenic_enrich <- fisher.test(intergenic_confusion, alternative = "greater")
    intergenic_enrich <- data.table( type="intergenic",p.value = ifelse(intergenic_enrich$p.value==0,2.2e-16, intergenic_enrich$p.value),  OR=intergenic_enrich$estimate)
    return(intergenic_enrich)
  }))
  message("      Intergenic hits: ",nrow(promoterHits))
  rm(txinfo)

  # combine all hits:
  snpHits <- do.call(rbind,list(cpgHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                enhancerHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                promoterHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                exonHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                # cdsHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                intronHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                utr3Hits[,c("chrom", "pos", "pValue", "anno", "type")],
                                utr5Hits[,c("chrom", "pos", "pValue", "anno", "type")],
                                tfHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                splicingHits[,c("chrom", "pos", "pValue", "anno", "type")],
                                intergenicHits[,c("chrom", "pos", "pValue", "anno", "type")]))
  # For the rest of the failed to be annotated, all in intergenic:
  failedSnps <- fsetdiff(snpInfo, snpHits[,.(chrom, pos, pValue)])
  failedSnps$anno <- "-"
  failedSnps$type <- "intergenic"
  snpHits <- rbind(snpHits, failedSnps)

  # prop:
  hitsProp <- round(prop.table(table(snpHits$type))*100, 3)
  hitsProp <- data.table(type=names(hitsProp), proportion=as.numeric(hitsProp))
  hitsProp <- hitsProp[order(-proportion)]
  message("==> Proportion: \n", paste("     ",paste(hitsProp$type, ":\t",hitsProp$proportion, "%",sep=""), collapse = "; \n"))

  # enrichment p-value and or:
  snpEnrich <- rbind(cpg_enrich, enhancer_enrich, promoter_enrich, exon_enrich, intron_enrich, utr3_enrich, utr5_enrich, tf_enrich, splicing_enrich,intergenic_enrich)

  return(list(snpHits = snpHits, snpEnrich = snpEnrich))
}

# snpHits <- xQTLanno_genomic(gwasDT[,.(chr, position,`p-value`)], p_cutoff = 5e-8, genomeVersion = "hg19")
# annoPlotBar <- xQTLvisual_anno(snpHits2,plotType="bar")






