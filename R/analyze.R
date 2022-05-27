#' @title xQTLanalyze_getSentinelSnp
#' @description detect sentinel SNPs in a given GWAS dataset.
#'  Return sentinel snps whose pValue < 5e-8(default) and SNP-to-SNP distance > 1e4 bp.
#' @param gwasDF A data.frame or a data.table object. Five columns are required (arbitrary column names is supported):
#'  `Col 1`. "snps" (character), , using an rsID (e.g. "rs11966562");
#'  `Col 2`. "chromosome" (character), one of the chromosome from chr1-chr22;
#'  `Col 3`. "postion" (integer), genome position of snp.
#'  `Col 4`. "P-value" (numeric).
#'  `Col 5`. "MAF" (numeric). Allel frequency.
#' @param pValueThreshold Cutoff of gwas p-value. Default: 5e-8
#' @param centerRange SNP-to-SNP distance. Default:1e6
#' @param mafThreshold Cutoff of maf to remove rare variants.
#' @param genomeVersion Genome version of input file. "grch37" or "grch38" (default).
#' @param grch37To38 TRUE or FALSE, we recommend converting grch37 to grch38, or using a input file of grch38 directly. Package `rtracklayer` is required.
#' @import data.table
#' @import stringr
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @importFrom GenomeInfoDb seqnames
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#'    gwasFile <- tempfile(pattern = "file")
#'    gwasURL <- paste0("https://raw.githubusercontent.com/dingruofan/",
#'                      "exampleData/master/gwas/AD/GLGC_AD_chr1_6_Sub3.txt")
#'    utils::download.file(gwasURL, destfile=gwasFile)
#'    gwasDF <- data.table::fread(gwasFile, sep="\t")
#'    gwasDF <- gwasDF[, .(rsid, chr, position, P, maf)]
#'    sentinelSnpDF <- xQTLanalyze_getSentinelSnp(gwasDF)
#' }
xQTLanalyze_getSentinelSnp <- function(gwasDF, pValueThreshold=5e-8, centerRange=1e6, mafThreshold = 0.01, genomeVersion="grch38", grch37To38 = FALSE){
  position <- pValue <- maf <- rsid <- chr <- NULL
  .<-NULL
  # Detect gene with sentinal SNP:
  gwasDF <- gwasDF[,1:5]
  message("== Preparing GWAS dataset... ")
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf")
  gwasDF <- na.omit(gwasDF)
  # chromosome revise：
  if( !str_detect(gwasDF[1,]$chr, regex("^chr")) ){
    gwasDF$chr <- paste0("chr", gwasDF$chr)
  }

  # convert variable class:
  gwasDF[,c("position", "pValue", "maf")] <- gwasDF[,.(position=as.numeric(position), pValue=as.numeric(pValue), maf=as.numeric(maf))]

  # MAF filter:
  gwasDF <- gwasDF[maf > mafThreshold & maf<1,]
  # 去重：
  # gwasDF <- gwasDF[order(rsid, pValue)][!duplicated(rsid)]
  # retain SNPs with rs id:
  # gwasDF <- gwasDF[stringr::str_detect(rsid,stringr::regex("^rs")),]

  ####################### convert hg19 to hg38:
  if(genomeVersion =="grch37" & grch37To38){
    message("== Converting SNPs' coordinate to GRCH38... ")

    if(suppressMessages(!requireNamespace("rtracklayer"))){
      message("Package [rtracklayer] is not installed! please install [rtracklayer] with following: ")
      message("---------")
      message('\"if (!require("BiocManager", quietly = TRUE)); BiocManager::install("rtracklayer")\"')
      message("---------")
    }

    path = system.file(package="xQTLbiolinks", "extdata", "hg19ToHg38.over.chain")
    ch = rtracklayer::import.chain(path)
    gwasRanges <- GenomicRanges::GRanges(gwasDF$chr,
                                         IRanges::IRanges(gwasDF$position, gwasDF$position),
                                         strand= "*",
                                         gwasDF[,.(rsid, maf, pValue)]
    )
    gwasRanges_hg38 <- unlist(rtracklayer::liftOver(gwasRanges, ch))
    # retain the SNPs located in chromosome 1:22
    gwasRanges_hg38 <- gwasRanges_hg38[which(as.character(GenomeInfoDb::seqnames(gwasRanges_hg38)) %in% paste0("chr", 1:22)),]
    gwasDF <- cbind(data.table(rsid = gwasRanges_hg38$rsid, maf=gwasRanges_hg38$maf, pValue=gwasRanges_hg38$pValue),
                    data.table::data.table(chr = as.character(GenomeInfoDb::seqnames(gwasRanges_hg38)), position=as.data.table(IRanges::ranges(gwasRanges_hg38))$start ))
    gwasDF <- gwasDF[,.(rsid, chr, position, pValue, maf)]
    message("== ",length(gwasRanges_hg38),"/",nrow(gwasDF)," left.")
    rm(gwasRanges, gwasRanges_hg38)
  }else if(genomeVersion =="grch37" & !grch37To38){
    stop("Please set \"grch37To38=TRUE\". Because the genome version of eqtl associations only support GRCH38, so sentinel SNP should be converted to GRCH38 before colocalization analysis if GRCH37 is provided.")
  }else if(genomeVersion =="grch38" & grch37To38){
    stop("Only grch37 genome version can be converted to grch38!")
  }else if(!genomeVersion %in% c("grch38", "grch37")){
    stop("Paramater genomeVersion must be choosen from grch38 and grch37, default: grch38.")
  }

  ####################### detect sentinel snp:
  message("== Detecting sentinel SNPs... ")
  # pValue filter:
  gwasDFsub <- gwasDF[pValue<pValueThreshold, ]
  chrAll <- unique(gwasDFsub$chr)
  chrAll <- chrAll[order(as.numeric(str_remove(chrAll,"chr")))]
  sentinelSnpDF <- data.table()
  for(i in 1:length(chrAll)){
    gwasDFsubChrom <- gwasDFsub[chr==chrAll[i],][order(pValue, position)]
    tmp <- copy(gwasDFsubChrom)
    sentinelSnps_count <- 0
    while(nrow(tmp)>0){
      sentinelSnps_count <- sentinelSnps_count+1
      sentinelSnpDF <- rbind(sentinelSnpDF, tmp[1,])
      startPos <- tmp[1,]$position-centerRange
      endPos <- tmp[1,]$position+centerRange
      tmp <- tmp[position<startPos | position>endPos][order(pValue, position)]
    }
    message("   In ",chrAll[i], ", ",sentinelSnps_count, " sentinel SNPs detected. ")
    rm(tmp, sentinelSnps_count, startPos, endPos)
  }
  message("== Totally, [", nrow(sentinelSnpDF), "] sentinel SNPs located in [",length(unique(sentinelSnpDF$chr))  ,"] chromosomes have been detected!")
  return(sentinelSnpDF)
}

#' @title xQTLanalyze_getTraits
#' @description identify trait genes with sentinel SNPs:
#' @param sentinelSnpDF A data.table. Better be the results from the function "xQTLanalyze_getSentinelSnp", five columns are required, including "rsid", "chr", "position", "pValue", and "maf".
#' @param detectRange A integer value. Trait genes that harbor sentinel SNPs located in the 1kb range upstream and downstream of gene. Default: 1e6 bp
#' @param tissueSiteDetail All tissues' name can be listed with \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param genomeVersion "grch38" or "grch37". Default: "grch38"
#' @param grch37To38 TRUE or FALSE, we recommend converting grch37 to grch38, or using a input file of grch38 directly. Package `rtracklayer` is required.
#' @import data.table
#' @import stringr
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @importFrom GenomeInfoDb seqnames
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#'   URL1<-"https://gitee.com/stronghoney/exampleData/raw/master/gwas/GLGC_CG0052/sentinelSnpDF.txt"
#'
#'   sentinelSnpDF <- data.table::fread(rawToChar(curl::curl_fetch_memory(URL1)$content))
#'   traitsAll <- xQTLanalyze_getTraits(sentinelSnpDF,detectRange=1e4,"Brain - Cerebellum",
#'                                      genomeVersion="grch37", grch37To38=TRUE)
#' }
xQTLanalyze_getTraits <- function(sentinelSnpDF, detectRange=1e6, tissueSiteDetail="", genomeVersion="grch38", grch37To38=FALSE){
  rsid <- maf <- strand <- pValue <- chr <- position <- chromosome <- NULL
  . <-genes <- geneSymbol <- gencodeId <- geneType <- description<- NULL

  data.table::as.data.table(sentinelSnpDF)


  if( length(tissueSiteDetail)!=1 | tissueSiteDetail=="" | !(tissueSiteDetail %in% tissueSiteDetailGTExv8$tissueSiteDetail) ){
    stop("== \"tissueSiteDetail\" can not be null. Please choose the tissue from tissue list of tissueSiteDetailGTExv8")
  }

  # (未做) 由于下一步的 xQTLdownload_eqtlPost 函数只能 query 基于 hg38(v26) 的突变 1e6 bp附近的基因，所以如果输入的GWAS是 hg19 的突变坐标，需要进行转换为38，然后再进行下一步 eqtl sentinel snp filter.
  # 由于从 EBI category 里获得的是 hg38(v26) 的信息，所以如果这一步是 hg19 的1e6范围内，则在 hg38里就会未必，所以需要这一步，如果是hg19，则对突变的坐标进行变换：
  ####################### convert hg19 to hg38:
  if( genomeVersion =="grch37" & !grch37To38 ){
    stop("Because the genome version of eqtl associations only support GRCH38, so sentinel SNP should be converted to GRCH38 before colocalization analysis if GRCH37 is provided.")
    genecodeVersion = "v19"
    datasetId="gtex_v7"
  }else if( genomeVersion =="grch37" & grch37To38){
    genecodeVersion = "v26"
    datasetId="gtex_v8"

    message("== Converting SNPs' coordinate to GRCH38... ")

    if(suppressMessages(!requireNamespace("rtracklayer"))){
      message("Package [rtracklayer] is not installed! please install [rtracklayer] with following: ")
      message("---------")
      message('\"if (!require("BiocManager", quietly = TRUE)); BiocManager::install("rtracklayer")\"')
      message("---------")
    }

    path = system.file(package="xQTLbiolinks", "extdata", "hg19ToHg38.over.chain")
    ch = rtracklayer::import.chain(path)
    dataRanges <- GenomicRanges::GRanges(sentinelSnpDF$chr,
                                         IRanges::IRanges(sentinelSnpDF$position, sentinelSnpDF$position),
                                         strand= "*",
                                         sentinelSnpDF[,.(rsid, maf, pValue)]
    )
    dataRanges_hg38 <- unlist(rtracklayer::liftOver(dataRanges, ch))
    # retain the SNPs located in chromosome 1:22
    dataRanges_hg38 <- dataRanges_hg38[which(as.character(GenomeInfoDb::seqnames(dataRanges_hg38)) %in% paste0("chr", 1:22)),]
    message("== ",length(dataRanges_hg38),"/",nrow(sentinelSnpDF)," left.")
    sentinelSnpDF <- cbind(data.table(rsid = dataRanges_hg38$rsid, maf=dataRanges_hg38$maf, pValue=dataRanges_hg38$pValue),
                    data.table::data.table(chr = as.character(GenomeInfoDb::seqnames(dataRanges_hg38)), position=as.data.table(IRanges::ranges(dataRanges_hg38))$start ))
    sentinelSnpDF <- sentinelSnpDF[,.(rsid, chr, position, pValue, maf)]
    rm(dataRanges, dataRanges_hg38)
  }else if( genomeVersion =="grch38" & !grch37To38){
    genecodeVersion = "v26"
    datasetId="gtex_v8"
  }else if( genomeVersion =="grch38" & grch37To38 ){
    genecodeVersion = "v26"
    datasetId="gtex_v8"
    message("Genome version of grch38 is used!")
  }else if(!genomeVersion %in% c("grch38", "grch37")){
    stop("Paramater genomeVersion must be choosen from grch38 and grch37, default: grch38.")
  }



  ###################### tissueSiteDetail="Brain - Cortex", maf_threshold=0.01
  geneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges, genomeVersion = genecodeVersion)
  chrAll <- unique(sentinelSnpDF$chr)
  chrAll <- chrAll[order(as.numeric(str_remove(chrAll,"chr")))]

  ########################## # 确保这些 SNPs 在 GTEx 里是有的， rsid 和 位置信息都一致才保留：
  # message("== Validating variant in GTEx....")
  # sentinelSnpDFNew <- data.table()
  # for(i in 1:length(chrAll)){
  #   sentinelSnpChromRaw <- sentinelSnpDF[chr==chrAll[i],]
  #   snpTMP <- suppressMessages(xQTLquery_varPos(chrom = ifelse(datasetId=="gtex_v8", chrAll[i], str_remove(chrAll[i], "chr")), pos = sentinelSnpChromRaw$position, datasetId = datasetId))
  #   if( !exists("snpTMP") || is.null(snpTMP)  ){
  #     stop("== Fail to fecth variant information, please check your network.")
  #   }
  #   # 使用 position 和 rsid 进行merge:， 因为 rsid 有可能有版本差异，而position可能有由于 liftover 转换产生差异。
  #   sentinelSnpChrom <- merge(sentinelSnpChromRaw, unique(snpTMP[,.(position=pos,rsid = snpId)]), by=c("rsid", "position"))
  #   if(nrow(sentinelSnpChrom)>0){
  #     sentinelSnpDFNew <- rbind(sentinelSnpDFNew, sentinelSnpChrom)
  #   }else{
  #     next()
  #   }
  #   message("   ", chrAll[i], ", ", nrow(sentinelSnpChrom),"/", nrow(sentinelSnpChromRaw) ," SNPs retained.")
  #   rm(sentinelSnpChromRaw, sentinelSnpChrom, snpTMP)
  # }
  # if(nrow(sentinelSnpDFNew)==0){
  #   stop("There is no shared variants found in your gwas dataset.")
  # }else{
  #   message("== Totally, ",nrow(sentinelSnpDFNew), "/",nrow(sentinelSnpDF), " sentinel SNPs retained.")
  # }

  sentinelSnpDFNew <- copy(sentinelSnpDF)

  ######################## 1
  message("== Detecting the trait genes of sentinel snps... ")
  traitsAll <- data.table()
  for(i in 1:length(chrAll)){
    geneInfoChrom <- geneInfo[chromosome == chrAll[i]]
    sentinelSnpChrom <- sentinelSnpDFNew[chr==chrAll[i],]
    Traits <- rbindlist(lapply(1:nrow(geneInfoChrom), function(x){
      tmp<-geneInfoChrom[x,];
      startPos <- tmp$start - detectRange; endPos <- tmp$end + detectRange;
      traits <- sentinelSnpChrom[position>=startPos & position<=endPos ];
      traits$gencodeId <- tmp$gencodeId
      return(traits)
    }))
    traitsAll <- rbind(traitsAll, Traits)
    message(i," - ", chrAll[i]," - [",nrow(Traits),"] gene-SNP pairs of [",length(unique(Traits$gencodeId)),"] genes and [",length(unique(Traits$rsid)),"] SNPs." )
    rm(sentinelSnpChrom, geneInfoChrom, Traits)
  }
  message("== Fetching [", length(unique(traitsAll$gencodeId)),"] genes' information from the GTEx.")
  geneAnnot <- xQTLquery_gene(unique(traitsAll$gencodeId), geneType = "gencodeId", gencodeVersion = genecodeVersion)
  if( datasetId=="gtex_v7" ){
    geneAnnot$chromosome <- paste0("chr",geneAnnot$chromosome)
  }
  if( exists("geneAnnot") && !is.null(geneAnnot) && nrow(geneAnnot)>0 ){
    traitsAll <- merge(geneAnnot[,.(genes, geneSymbol,gencodeId, geneType, description, chromosome, start,end ,strand)],traitsAll, by.x="genes",by.y="gencodeId",  all.x=TRUE)[,-c("genes")]
    traitsAll <- traitsAll[,.( chromosome, geneStart=start, geneEnd=end, geneStrand=strand, geneSymbol, gencodeId, rsid, position, pValue, maf )][order(as.numeric(str_remove(chromosome, "chr")), pValue, position)]
    message("== Totally, [",nrow(traitsAll), "] associations between [",length(unique(traitsAll$gencodeId)),"] traits genes and [",length(unique(traitsAll$rsid)),"] SNPs are detected." )
  }

  # Get the overlap with the eGgenes:
  egeneDF <- xQTLdownload_egene(tissueSiteDetail = tissueSiteDetail) #11240
  traitsAll <- traitsAll[gencodeId %in% egeneDF$gencodeId]
  message("== After taking the intersection with egenes, [",nrow(traitsAll), "] associations between [",length(unique(traitsAll$gencodeId)),"] traits genes and [",length(unique(traitsAll$rsid)),"] SNPs are detected." )

  return(traitsAll)
}


#' @title xQTLanalyze_coloc
#' @description conduct colocalization analysis with detected gene.
#'
#' @param gwasDF A dataframe of gwas.
#' @param traitGene A gene symbol or a gencode id (versioned).
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param genomeVersion "v26" (default) or "v19"
#' @param tissueSiteDetail A character string. Tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param mafThreshold Cutoff of maf to remove rare variants.
#' @param population Supported population is consistent with the LDlink, which can be listed using function "LDlinkR::list_pop()"
#' @param gwasSampleNum Sample number of GWAS dataset. Default:50000.
#' @param token LDlink provided user token, default = NULL, register for token at https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param method Now only one "coloc". Package `coloc` is required.
#'
#' @return coloc resut
#' @export
#'
#' @examples
#' \donttest{
#'   # please see see vignette.
#' }
xQTLanalyze_coloc <- function(gwasDF, traitGene, geneType="auto", genomeVersion="grch38", tissueSiteDetail="", mafThreshold=0.01, population="EUR", gwasSampleNum=50000, method="coloc", token="9246d2db7917"){
  rsid <- chr <- position <- se <- pValue <- snpId <- maf <- i <- variantId <- NULL
  . <- NULL

  # tissueSiteDetail="Brain - Cortex"
  # geneType="geneSymbol"
  # genomeVersion="grch37"
  # gwasSampleNum=50000
  # mafThreshold=0.01
  population <- ""
  token <- ""

  if(!requireNamespace("coloc")){
    stop("please install package \"coloc\" with install.packages(\"coloc\").")
  }

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(traitGene, function(g){ stringr::str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  ###################### eqtl dataset:
  if( genomeVersion=="grch37"){
    gencodeVersion <- "v19"
    geneInfo <- xQTLquery_gene(traitGene, geneType = geneType, gencodeVersion =gencodeVersion)[1,]
    eqtlInfo <- xQTLdownload_eqtlAllAsso(traitGene,geneType = geneType, tissueSiteDetail=tissueSiteDetail, withB37VariantId = FALSE)
    # eqtlInfo <- eqtlInfo[b37VariantId!=""]
    eqtlInfo[,"position":= .( as.integer(unlist(lapply(variantId, function(x){str_split(x, fixed("_"))[[1]][2]}))) )]
    # chromosome:
    P_chrom <- paste0("chr",geneInfo$chromosome)
  }else if(genomeVersion=="grch38"){
    gencodeVersion <- "v26"
    geneInfo <- xQTLquery_gene(traitGene, geneType = geneType, gencodeVersion =gencodeVersion)[1,]
    eqtlInfo <- xQTLdownload_eqtlAllAsso(traitGene,geneType = geneType, tissueSiteDetail=tissueSiteDetail, withB37VariantId = FALSE)
    eqtlInfo[,"position":= .( unlist(lapply(variantId, function(x){str_split(x, fixed("_"))[[1]][2]})) )]
    # chromosome:
    P_chrom <- geneInfo$chromosome
  }

  if( !exists("eqtlInfo") || is.null(eqtlInfo) || nrow(eqtlInfo)==0){
    message(i," | gene", traitGene, "has no eqtl associations, next!")
    message(" = None eQTL associations obtained of gene [",traitGene,"], please change the gene name or ENSEMBLE ID.")
    return(NULL)
  }


  eqtlInfo<- eqtlInfo[maf>mafThreshold & maf <1,]
  # 去重：
  eqtlInfo <- eqtlInfo[order(snpId, pValue)][!duplicated(snpId)]
  # subset:
  eqtlInfo <- na.omit( eqtlInfo[,.(rsid=snpId, maf=as.numeric(maf), beta, se, pValue=as.numeric(pValue), position=as.numeric(position))] )

  message("Data processing...")

  ##################### gwas dataset:
  gwasDF <- gwasDF[,1:5]
  data.table::setDT(gwasDF)
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf")
  # gwas subset:
  gwasDF <- na.omit(gwasDF)

  gwasDF <- gwasDF[chr==ifelse(stringr::str_detect(gwasDF[1,]$chr, "chr"),P_chrom, stringr::str_remove(P_chrom, "chr")),]
  # convert variable class:
  gwasDF[,c("position", "pValue", "maf")] <- gwasDF[,.(position=as.numeric(position), pValue=as.numeric(pValue), maf=as.numeric(maf))]
  # MAF filter:
  gwasDF <- gwasDF[maf > mafThreshold & maf<1,]
  # 去重：
  gwasDF <- gwasDF[order(rsid, pValue)][!duplicated(rsid)]
  # retain SNPs with rs id:
  gwasDF <- gwasDF[stringr::str_detect(rsid,stringr::regex("^rs")),]
  # chromosome revise：
  if( !str_detect(gwasDF[1,]$chr, stringr::regex("^chr")) ){
    gwasDF$chr <- paste0("chr", gwasDF$chr)
  }

  ##### convert to grch38:
  if(genomeVersion == "grch37"){
    message("== Converting GWAS coordinate to GRCH38... ")

    if(suppressMessages(!requireNamespace("rtracklayer"))){
      message("Package [rtracklayer] is not installed! please install [rtracklayer] with following: ")
      message("---------")
      message('\"if (!require("BiocManager", quietly = TRUE)); BiocManager::install("rtracklayer")\"')
      message("---------")
    }

    path = system.file(package="xQTLbiolinks", "extdata", "hg19ToHg38.over.chain")
    ch = rtracklayer::import.chain(path)
    dataRanges <- GenomicRanges::GRanges(gwasDF$chr,
                                         IRanges::IRanges(gwasDF$position, gwasDF$position),
                                         strand= "*",
                                         gwasDF[,.(rsid, maf, pValue)]
    )
    dataRanges_hg38 <- unlist(rtracklayer::liftOver(dataRanges, ch))
    # retain the SNPs located in chromosome 1:22
    dataRanges_hg38 <- dataRanges_hg38[which(as.character(GenomeInfoDb::seqnames(dataRanges_hg38)) %in% paste0("chr", 1:22)),]
    message("== ",length(dataRanges_hg38),"/",nrow(gwasDF)," (",round(length(dataRanges_hg38)/nrow(gwasDF)*100,2),"%)"," left.")
    gwasDFnew <- cbind(data.table(rsid = dataRanges_hg38$rsid, maf=dataRanges_hg38$maf, pValue=dataRanges_hg38$pValue),
                       data.table::data.table(chr = as.character(GenomeInfoDb::seqnames(dataRanges_hg38)), position=as.data.table(IRanges::ranges(dataRanges_hg38))$start ))
    gwasDF <- gwasDFnew[chr==P_chrom,.(rsid, chr, position, pValue, maf)]
    rm(dataRanges, dataRanges_hg38, gwasDFnew)
  }


  #
  tissueSiteDetailId <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
  gwasEqtlInfo <- merge(gwasDF, eqtlInfo[,.(rsid, maf, pValue, position)], by=c("rsid", "position"), suffixes = c(".gwas",".eqtl"))
  # centerSnp <- gwasEqtldata[which.min(gwasEqtldata$pValue.gwas),]
  # centerSnp <- xQTLquery_varId(sentinelSnp, variantType = variantType)
  # if(nrow(centerSnp)==0 ){
  #   centerSnp <- gwasDF[which.min(gwasDF$pValue),]
  #   gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$pos-(colocRange)) & position<=(centerSnp$pos+(colocRange)),][order(position)]
  # }
  # if( colocRange==0 ){
  #   gwasEqtlInfo <- gwasEqtldata[order(position)]
  # }else {
  #   gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$pos-(colocRange)) & position<=(centerSnp$pos+(colocRange)),][order(position)]
  # }
  # if(nrow(gwasEqtlInfo)==0){
  #   centerSnp <- gwasDF[which.min(gwasDF$pValue),]
  #   gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$pos-(colocRange)) & position<=(centerSnp$pos+(colocRange)),][order(position)]
  # }
  if(nrow(gwasEqtlInfo)==0){
    message("No shared variants between eQTL and GWAS, please check your input!.")
    return(NULL)
  }
  message("== Start the colocalization analysis of gene [", traitGene,"]")
  if(suppressMessages(!requireNamespace("coloc"))){
    message("Package [coloc] is not installed! please install [coloc] with following: ")
    message("---------")
    message('\"if (!require("BiocManager", quietly = TRUE)); install.packages(\"coloc\")')
    message("---------")
  }
  # 只选择中心附近的 Num snp 进行分析：
  # 防止 check_dataset中 p = pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2 出错
  suppressWarnings(coloc_Out <- coloc::coloc.abf(dataset1 = list( pvalues = gwasEqtlInfo$pValue.gwas, type="quant", N=gwasSampleNum, snp=gwasEqtlInfo$rsid, MAF=gwasEqtlInfo$maf.gwas),
                                                 dataset2 = list( pvalues = gwasEqtlInfo$pValue.eqtl, type="quant", N=sampleNum[tissueSiteDetailId, on="tissueSiteDetailId"]$sampleNum, snp=gwasEqtlInfo$rsid, MAF= gwasEqtlInfo$maf.eqtl)))
  coloc_Out_results <- as.data.table(coloc_Out$results)
  # coloc_Out_results$gene <- traitGenes[i]
  coloc_Out_summary <- as.data.table(t(as.data.frame(coloc_Out$summary)))
  coloc_Out_summary$traitGene <- traitGene
  message("== Done")

  print(coloc_Out_summary)
  # coloc_Out_summary$pearsonCoor <- cor(-log(gwasEqtlInfo$pValue.gwas, 10),-log(gwasEqtlInfo$pValue.eqtl, 10), method = "pearson")

  return(list(coloc_Out_summary=coloc_Out_summary,coloc_Out_results=coloc_Out_results, gwasEqtlInfo=gwasEqtlInfo))
}


#' @title xQTLanalyze_TSExp
#' @description perform tissue-specific expression analysis.
#'
#' @param genes A charater vector or a string of gene symbol, gencode id (versioned), or a charater string of gene type.
#' @param geneType A character string. "auto"(default), "geneSymbol", "gencodeId" or "geneCategory".
#' @param method "SPM" or "entropy"
#' @param datasetId "gtex_v8" or "gtex_v7".
#'
#' @return A data.table
#' @export
#'
#' @examples
#' \donttest{
#'  TSgene <- xQTLanalyze_TSExp(extractGeneInfo(gencodeGeneInfoAllGranges)$gencodeId[1:20])
#'  # xQTLvisual_geneExpTissues( TSgene[order(-DPM)][1,]$geneSymbol )
#' }
xQTLanalyze_TSExp <- function(genes, geneType="auto", method="SPM", datasetId="gtex_v8"){
  gencodeId <- geneSymbol <-.<-NULL
  if(datasetId == "gtex_v8"){
    genomeVersion="v26"
    tissueSiteDetail <- copy(tissueSiteDetailGTExv8)
  }else if(datasetId == "gtex_v7"){
    genomeVersion="v19"
    tissueSiteDetail <- copy(tissueSiteDetailGTExv7)
  }
  geneExp <- xQTLdownload_geneMedExp(genes=genes,  geneType=geneType, datasetId = datasetId)
  geneExpCast <- dcast(geneExp[,.(gencodeId, geneSymbol,geneType, median, tissueSiteDetail)], gencodeId+geneSymbol+geneType~tissueSiteDetail, value.var = "median")
  geneExpCastMat <- geneExpCast[, -c("gencodeId", "geneSymbol", "geneType")]

  if( method=="SPM"){
    DPMlist <- lapply( 1:nrow(geneExpCastMat), function(x){
      x <- unlist(geneExpCastMat[x,])
      x<- as.numeric(x)
      if(all(x==0)){
        return(list(SPM=x, DPM=NA ))
      }else{
        SPMi <- unlist(lapply(1:length(x), function(xx){
          xi <- rep(0, length(x))
          xi[xx] <- x[xx]
          if(all(xi==0)){
            cosX = 0
          }else{
            # cosX <- xi%*%x/( length(xi)*length(x) )
            cosX <- crossprod(x, xi)/sqrt(crossprod(x) * crossprod(xi))
          }
          return(cosX)
        }))
        return(list(SPM=SPMi, DPM=sd(SPMi)*sqrt(length(x)) ))
      }
    })
    SPMmat <- as.data.table(t(as.data.table( lapply(DPMlist, function(x){ x[[1]] }))))
    names(SPMmat) <- names(geneExpCastMat)
    DPM <- unlist(lapply(DPMlist, function(x){ x[[2]] }))
    SPMmat$DPM <- DPM
    SPMmat <- cbind(geneExpCast[, c("gencodeId", "geneSymbol", "geneType")], SPMmat)

    return(SPMmat)
  }else if(method== "entropy"){
    TSI <- copy(geneExpCastMat)
    TSI$TSI <- apply(geneExpCastMat, 1, function(x){
      x <- as.numeric(x)
      if( all(x==0)){
        pi <- rep(0, length(x))
      }else{
        pi <- x/sum(x)
      }
      sum(na.omit(-pi*log(pi,2)))
    })
    TSI <- cbind(geneExpCast[, c("gencodeId", "geneSymbol", "geneType")], TSI)
    return(TSI)
  }else{
    stop("Please choose the right method.")
  }

}












