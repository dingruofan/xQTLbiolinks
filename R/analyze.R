#' @title Detect sentinel SNPs in a given summary statistis dataset.
#' @description
#'  Return sentinel snps whose pValue < 5e-8(default) and SNP-to-SNP distance > 1e6 bp.
#' @param gwasDF A data.frame or a data.table object. Five columns are required (arbitrary column names is supported):
#'
#'  `Col 1`. "snps" (character), , using an rsID (e.g. "rs11966562").
#'
#'  `Col 2`. "chromosome" (character), one of the chromosome from chr1-chr22.
#'
#'  `Col 3`. "postion" (integer), genome position of snp.
#'
#'  `Col 4`. "P-value" (numeric).
#'
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
#' url<-"http://raw.githubusercontent.com/dingruofan/exampleData/master/GLGC.txt"
#' gwasDF <- data.table::fread(url)
#' gwasDF <- gwasDF[, .(rsid, chr, position, P, maf)]
#' sentinelSnpDF <- xQTLanalyze_getSentinelSnp(gwasDF)
xQTLanalyze_getSentinelSnp <- function(gwasDF, pValueThreshold=5e-8, centerRange=1e6, mafThreshold = 0.01, genomeVersion="grch38", grch37To38 = FALSE){
  position <- pValue <- maf <- rsid <- chr <- NULL
  .<-NULL
  # Detect gene with sentinal SNP:
  gwasDF <- gwasDF[,1:5]
  data.table::setDT(gwasDF)
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
    stop("Please set \"grch37To38=TRUE\". Because the genome version of eqtl associations only support GRCH38, sentinel SNP should be converted to GRCH38 before colocalization analysis if GRCH37 is provided.")
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

#' @title Identify trait genes using sentinel SNPs generated from `xQTLanalyze_getSentinelSnp`
#' @description
#'
#' @param sentinelSnpDF A data.table. Better be the results from the function "xQTLanalyze_getSentinelSnp", five columns are required, including "rsid", "chr", "position", "pValue", and "maf".
#' @param detectRange A integer value. Trait genes that harbor sentinel SNPs located in the 1kb range upstream and downstream of gene. Default: 1e6 bp
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using "tissueSiteDetailGTExv8" or "tissueSiteDetailGTExv7"
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
#' URL1<-"http://gitee.com/stronghoney/exampleData/raw/master/gwas/GLGC_CG0052/sentinelSnpDF.txt"
#' sentinelSnpDF <- data.table::fread(URL1)
#' traitsAll <- xQTLanalyze_getTraits(sentinelSnpDF,detectRange=1e4,"Brain - Cerebellum",
#'                                    genomeVersion="grch37", grch37To38=TRUE)
xQTLanalyze_getTraits <- function(sentinelSnpDF, detectRange=1e6, tissueSiteDetail="", genomeVersion="grch38", grch37To38=FALSE){
  rsid <- maf <- strand <- pValue <- chr <- position <- chromosome <- NULL
  . <-genes <- geneSymbol <- gencodeId <- geneType <- description<- NULL

  data.table::setDT(sentinelSnpDF)


  if( length(tissueSiteDetail)!=1 | tissueSiteDetail=="" | !(tissueSiteDetail %in% tissueSiteDetailGTExv8$tissueSiteDetail) ){
    stop("== \"tissueSiteDetail\" can not be null. Please choose the tissue from tissue list of tissueSiteDetailGTExv8")
  }

  # (未做) 由于下一步的 xQTLdownload_eqtlPost 函数只能 query 基于 hg38(v26) 的突变 1e6 bp附近的基因，所以如果输入的GWAS是 hg19 的突变坐标，需要进行转换为38，然后再进行下一步 eqtl sentinel snp filter.
  # 由于从 EBI category 里获得的是 hg38(v26) 的信息，所以如果这一步是 hg19 的1e6范围内，则在 hg38里就会未必，所以需要这一步，如果是hg19，则对突变的坐标进行变换：
  ####################### convert hg19 to hg38:
  if( genomeVersion =="grch37" & !grch37To38 ){
    stop("Because the genome version of eqtl associations only support GRCH38, sentinel SNP should be converted to GRCH38 before colocalization analysis if GRCH37 is provided.")
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


#' @title Conduct colocalization analysis with trait genes generated from `xQTLanalyze_getTraits`
#' @description
#'
#' @param gwasDF A dataframe of gwas.
#' @param traitGene A gene symbol or a gencode id (versioned).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param genomeVersion "grch38" (default) or "grch37". Note: grch37 will be converted to grch38 automatically.
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using "tissueSiteDetailGTExv8" or "tissueSiteDetailGTExv7"
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
#' # Please see see vignette: https://dingruofan.github.io/xQTLbiolinks/articles/
xQTLanalyze_coloc <- function(gwasDF, traitGene, geneType="auto", genomeVersion="grch38", tissueSiteDetail="", mafThreshold=0.01, population="EUR", gwasSampleNum=50000, method="coloc", token="9246d2db7917"){
  rsid <- chr <- position <- se <- pValue <- snpId <- maf <- pos <- i <- variantId <- NULL
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
  ###################### eqtl dataset:
  eqtlInfo <- xQTLdownload_eqtlAllAsso(traitGene, geneType = geneType, tissueSiteDetail=tissueSiteDetail, withB37VariantId = FALSE)
  eqtlInfo[,position:=.(pos)]

  if( !exists("eqtlInfo") || is.null(eqtlInfo)){
    message(i," | gene", traitGene, "has no eqtl associations, next!")
    message(" = None eQTL associations obtained of gene [",traitGene,"], please change the gene name or ENSEMBLE ID.")
    return(NULL)
  }

  # chromosome:
  P_chrom <- paste0("chr",eqtlInfo[1,]$chrom)

  eqtlInfo<- eqtlInfo[maf>mafThreshold & maf <1,]
  # 去重：
  eqtlInfo <- eqtlInfo[order(snpId, pValue)][!duplicated(snpId)]
  # subset:
  eqtlInfo <- na.omit( eqtlInfo[,.(rsid=snpId, maf=as.numeric(maf), beta, se, pValue=as.numeric(pValue), position=as.numeric(position))] )

  if(nrow(eqtlInfo)==0){
    message("Number of eQTL association is <1")
    return(NULL)
  }



  message("Data processing...")

  ##################### gwas dataset:
  gwasDF <- gwasDF[,1:5]
  data.table::setDT(gwasDF)
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf")
  # gwas subset:
  gwasDF <- na.omit(gwasDF)

  p_chrom_tmp <- ifelse(stringr::str_detect(gwasDF[1,]$chr, "chr"),P_chrom, stringr::str_remove(P_chrom, "chr"))
  gwasDF <- gwasDF[chr==p_chrom_tmp,]
  # convert variable class:
  gwasDF[,c("position", "pValue", "maf"):=.(as.numeric(position), as.numeric(pValue), as.numeric(maf))]
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


#' @title Perform tissue-specific expression analysis for genes.
#' @description
#'
#' @param genes A charater vector or a string of gene symbol, gencode id (versioned), or a charater string of gene type.
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param method "SPM" or "entropy"
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#'
#' @return A data.table
#' @export
#'
#' @examples
#' TSgene <- xQTLanalyze_TSExp(extractGeneInfo(gencodeGeneInfoAllGranges)$gencodeId[1:20])
#' # xQTLvisual_geneExpTissues( TSgene[order(-DPM)][1,]$geneSymbol )
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





#' @title eQTL-specific analysis
#'
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param binNum A integer value. Number of bins to split values of R2 of LD into homogeneous bins.
#' @param study Studies can be listed using "ebi_study_tissues"
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' p <- xQTLanalyze_qtlSpecificity(gene="MMP7", variantName="rs11568818", study="gtex_v8")
xQTLanalyze_qtlSpecificity <- function(gene="", geneType="auto", variantName="", variantType="auto", binNum=4, study="", population="EUR"){

  ebi_ST <-copy(ebi_study_tissues)
  if(gene=="" || variantName==""){
    stop("gene and variant can not be null!")
  }

  # check geneType
  if( !(geneType %in% c("auto","geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"auto\", \"geneSymbol\", and \"gencodeId\".")
  }
  if( length(gene)==1 && gene!=""){
    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }
  }

  # check variantType:
  if( !(variantType %in% c("auto","snpId", "variantId")) ){
    stop("Parameter \"geneType\" should be choosen from \"auto\", \"snpId\", and \"variantId\".")
  }
  if(length(variantName)==1 && variantName!=""){
    # auto pick variantType
    if(variantType=="auto"){
      if(stringr::str_detect(variantName, stringr::regex("^rs"))){
        variantType <- "snpId"
      }else if(stringr::str_count(variantName,"_")>=3){
        variantType <- "variantId"
      }else{
        stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
      }
    }
  }

  # check study:
  if( length(study) ==1 && study!="" ){
    if(toupper(study) %in% toupper(unique(ebi_ST$study_accession))){
      study <- unique(ebi_ST$study_accession)[ toupper(unique(ebi_ST$study_accession)) == toupper(study) ]
      message("== Study [", study, "] detected...")
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== Study [",study,"] can not be correctly matched, please choose from above list: ")
    }
  }

  # variant detail:
  variantInfo <- xQTLquery_varId(variantName, variantType = variantType)

  # fetch LD:
  message("== Retrieve LD information of SNP: [",variantInfo$snpId,"]...")
  try(snpLD <- retrieveLD(variantInfo$chromosome, variantInfo$snpId, population))
  # try(snpLD <- retrieveLD_LDproxy(highlightSnp,population = population, windowSize = windowSize, genomeVersion = genomeVersion, token = token) )
  data.table::setDT(snpLD)

  if(nrow(snpLD)<1){
    stop("No LD found!")
  }

  # cut LD into bins:
  snpLD$LDbins <- as.character(cut(snpLD$R2, breaks=seq(0,1,length.out=(10+1)) ))
  snpLD <- snpLD[order(R2)]

  message("== Start download associations of QTL...")
  assoAll <- data.table()
  for(i in 1:nrow(snpLD)){
    suppressMessages( asso_I <- xQTLdownload_eqtlAllAsso(gene = gene, variantName= snpLD[i,]$SNP_B, study=study, withB37VariantId=FALSE) )
    if( is.null(asso_I) || nrow(asso_I)==0){
      message("== Num:",i, ", SNP: [",snpLD[i,]$SNP_B,"], got records: ", 0,", Skipped. ", format(Sys.time(), "| %Y-%b-%d %H:%M:%S "))
      next()
    }else{
      message("== Num:",i, "; SNP: [",snpLD[i,]$SNP_B,"]; got records: [", nrow(asso_I),"]; Tissues: [", length(unique(asso_I$tissue_label)),"]; Studies: [", length(unique(asso_I$study_id)), format(Sys.time(), "] | %Y-%b-%d %H:%M:%S "))
    }
    asso_I <- asso_I[,.(snpId, pValue, tissue, tissue_label, study_id, qtl_group)]
    assoAll <- rbind(assoAll,asso_I)
    rm(asso_I)
  }
  # Retain min pvalue in each tissue_label-study, due to the duplication induced by qtl_group.
  assoAll <- assoAll[,.SD[which.min(pValue),], by=c("tissue", "tissue_label", "study_id", "qtl_group", "snpId")]
  assoAllLd <- merge(snpLD[,.(snpId=SNP_B, LDbins, R2)], assoAll, by="snpId")
  assoAllLd$logP <- ifelse(assoAllLd$pValue==0, 0, (-log(assoAllLd$pValue,10)))
  # Retain min/max pvalue in each bin:
  minP_f <- function(x){ x[logP==max(x$logP),][1,] }
  assoAllLd <- assoAllLd[,minP_f(.SD), by=c("LDbins", "tissue", "tissue_label", "study_id", "qtl_group")]
  # assoAllLd <- merge(assoAllLd, ebi_ST[,.(tissueLabel=tissue_label, tissue)], by="tissue")

  assoAllLd <- merge(assoAllLd, assoAllLd[,.(corRP=cor(R2, logP)),by="tissue_label"][order(corRP)], by="tissue_label")

  ggplot(assoAllLd)+
    geom_tile(aes(x=LDbins, y=reorder(tissue_label, corRP), fill=logP),color = "black")+
    geom_text(aes(x=LDbins, y=tissue_label, label = round(logP,2),  color = logP), size = 3.5)+
    scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
    theme_classic()

  plot(p)
  return(p)

}










