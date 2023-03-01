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
#'
#'  `Col 6`. "beta" (numeric). effect size.
#'
#'  `Col 7`. "se" (numeric). standard error.
#'
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
#' url<-"http://raw.githubusercontent.com/dingruofan/exampleData/master/GLGC.txt"
#' gwasDF <- data.table::fread(url)
#' gwasDF <- gwasDF[, .(rsid, chr, position, P, maf, beta, se)]
#' sentinelSnpDF <- xQTLanalyze_getSentinelSnp(gwasDF)
#' }
xQTLanalyze_getSentinelSnp <- function(gwasDF, pValueThreshold=5e-8, centerRange=1e6, mafThreshold = 0.01, genomeVersion="grch38", grch37To38 = FALSE){
  position <- pValue <- maf <- rsid <- chr <- se <- NULL
  .<-NULL
  # Detect gene with sentinal SNP:
  gwasDF <- gwasDF[,1:7]
  data.table::setDT(gwasDF)
  message("== Preparing GWAS dataset... ")
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf", "beta", "se")
  gwasDF <- na.omit(gwasDF)
  # chromosome revise：
  if( !str_detect(gwasDF[1,]$chr, regex("^chr")) ){
    gwasDF$chr <- paste0("chr", gwasDF$chr)
  }

  # convert variable class:
  gwasDF[,c("position", "pValue", "maf", "beta", "se")] <- gwasDF[,.(position=as.numeric(position), pValue=as.numeric(pValue), maf=as.numeric(maf), beta=as.numeric(beta), se=as.numeric(se))]

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
                                         gwasDF[,.(rsid, maf, pValue, beta, se)]
    )
    gwasRanges_hg38 <- unlist(rtracklayer::liftOver(gwasRanges, ch))
    # retain the SNPs located in chromosome 1:22
    gwasRanges_hg38 <- gwasRanges_hg38[which(as.character(GenomeInfoDb::seqnames(gwasRanges_hg38)) %in% paste0("chr", 1:22)),]
    gwasDF <- cbind(data.table(rsid = gwasRanges_hg38$rsid, maf=gwasRanges_hg38$maf, pValue=gwasRanges_hg38$pValue, beta = gwasRanges_hg38$beta, se = gwasRanges_hg38$se),
                    data.table::data.table(chr = as.character(GenomeInfoDb::seqnames(gwasRanges_hg38)), position=as.data.table(IRanges::ranges(gwasRanges_hg38))$start ))
    gwasDF <- gwasDF[,.(rsid, chr, position, pValue, maf, beta, se)]
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
#' @param sentinelSnpDF A data.table. Better be the results from the function "xQTLanalyze_getSentinelSnp", seven columns are required, including "rsid", "chr", "position", "pValue", "maf", "beta" and "se".
#' @param detectRange A integer value. Trait genes that harbor sentinel SNPs located in the 1kb range upstream and downstream of gene. Default: 1e6 bp
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param genomeVersion "grch38" or "grch37". Default: "grch38"
#' @param grch37To38 TRUE or FALSE, we recommend converting grch37 to grch38, or using a input file of grch38 directly. Package `rtracklayer` is required.
#' @param overlapWithEGene take the intersection with eGenes. Default:TRUE
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
#' URL1<-"https://gitee.com/stronghoney/exampleData/raw/master/gwas/GLGC_CG0052/sentinelSnpDF.txt"
#' sentinelSnpDF <- data.table::fread(URL1)
#' traitsAll <- xQTLanalyze_getTraits(sentinelSnpDF,detectRange=1e4,"Brain - Cerebellum",
#'                                    genomeVersion="grch37", grch37To38=TRUE)
#' }
xQTLanalyze_getTraits <- function(sentinelSnpDF, detectRange=1e6, tissueSiteDetail="", genomeVersion="grch38", grch37To38=FALSE, overlapWithEGene=TRUE){
  rsid <- maf <- strand <- pValue <- chr <- position <- chromosome <- NULL
  . <-genes <- geneSymbol <- gencodeId <- geneType <- description<- NULL

  data.table::setDT(sentinelSnpDF)

  if(overlapWithEGene){
    if( length(tissueSiteDetail)!=1 | tissueSiteDetail=="" | !(tissueSiteDetail %in% tissueSiteDetailGTExv8$tissueSiteDetail) ){
      stop("== \"tissueSiteDetail\" can not be null. Please choose the tissue from tissue list of tissueSiteDetailGTExv8")
    }
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
  geneAnnot <- xQTLquery_gene(unique(traitsAll$gencodeId), geneType = "gencodeId")
  if( datasetId=="gtex_v7" ){
    geneAnnot$chromosome <- paste0("chr",geneAnnot$chromosome)
  }
  if( exists("geneAnnot") && !is.null(geneAnnot) && nrow(geneAnnot)>0 ){
    traitsAll <- merge(geneAnnot[,.(genes, geneSymbol,gencodeId, geneType, description, chromosome, start,end ,strand)],traitsAll, by.x="genes",by.y="gencodeId",  all.x=TRUE)[,-c("genes")]
    traitsAll <- traitsAll[,.( chromosome, geneStart=start, geneEnd=end, geneStrand=strand, geneSymbol, gencodeId, rsid, position, pValue, maf )][order(as.numeric(str_remove(chromosome, "chr")), pValue, position)]
    message("== Totally, [",nrow(traitsAll), "] associations between [",length(unique(traitsAll$gencodeId)),"] traits genes and [",length(unique(traitsAll$rsid)),"] SNPs are detected." )
  }

  # Get the overlap with the eGgenes:
  if(overlapWithEGene){
    egeneDF <- xQTLdownload_egene(tissueSiteDetail = tissueSiteDetail) #11240
    traitsAll <- traitsAll[gencodeId %in% egeneDF$gencodeId]
    message("== After taking the intersection with egenes, [",nrow(traitsAll), "] associations between [",length(unique(traitsAll$gencodeId)),"] traits genes and [",length(unique(traitsAll$rsid)),"] SNPs are detected." )
  }

  return(traitsAll)
}


#' @title Conduct colocalization analysis with trait genes generated from `xQTLanalyze_getTraits`
#' @param gwasDF A data.frame or data.table objectof gwas.
#'
#' @param traitGene A gene symbol or a gencode id (versioned).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param genomeVersion "grch38" (default) or "grch37". Note: grch37 will be converted to grch38 automatically.
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param study (character) name of studies can be listed using "ebi_study_tissues"
#' @param mafThreshold Cutoff of maf to remove rare variants.
#' @param population Supported population is consistent with the LDlink, which can be listed using function "LDlinkR::list_pop()"
#' @param gwasSampleNum Sample number of GWAS dataset. Default:50000.
#' @param token LDlink provided user token, default = NULL, register for token at https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param method (character) options: "coloc"(default) or "hyprcoloc". Package `coloc` or `hyprcoloc` is required.
#' @param bb.alg For `hyprcoloc`, branch and bound algorithm: TRUE, employ BB algorithm; FALSE, do not. Default: FALSE.
#'
#' @return A list of coloc result and details.
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://raw.githubusercontent.com/dingruofan/exampleData/master/gwasDFsub_MMP7.txt"
#' gwasDF <- data.table::fread(url1)
#' output <- xQTLanalyze_coloc(gwasDF = gwasDF, traitGene= "MMP7", tissueSiteDetail="Prostate")
#' }
xQTLanalyze_coloc <- function(gwasDF, traitGene, geneType="auto", genomeVersion="grch38", tissueSiteDetail="", study="gtex_v8", mafThreshold=0.01, population="EUR", gwasSampleNum=50000, method="coloc", token="9246d2db7917", bb.alg=FALSE){
  rsid <- chr <- position <- se <- pValue <- snpId <- maf <- pos <- i <- variantId <- se.eqtl <- se.gwas <-SNP.PP.H4 <- beta.eqtl <- beta.gwas <- posterior_prob <- regional_prob <-candidate_snp <- posterior_explained_by_snp  <- NULL
  qtl_group <- NULL
  . <- NULL

  # tissueSiteDetail="Brain - Cortex"
  # geneType="geneSymbol"
  # genomeVersion="grch37"
  # gwasSampleNum=50000
  # mafThreshold=0.01
  # https://stackoverflow.com/questions/66849936/make-cran-r-package-suggest-github-r-package
  population <- ""
  token <- ""

  if(method == "coloc" && !requireNamespace("coloc")){
    stop("please install package \"coloc\" with install.packages(\"coloc\").")
  }

  if(method == "hyprcoloc" && !requireNamespace("hyprcoloc")){
    stop("please install package \"hyprcoloc\".")
  }

  if(method == "Both"){
    if( !requireNamespace("coloc") ){ stop("please install package \"coloc\" with install.packages(\"coloc\").") }
    if( !requireNamespace("hyprcoloc") ){ stop("please install package \"hyprcoloc\".") }
  }

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(traitGene, function(g){ stringr::str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  # eqtl dataset:
  eqtlInfo <- xQTLdownload_eqtlAllAsso(traitGene, geneType = geneType, tissueLabel=tissueSiteDetail, study=study, withB37VariantId = FALSE)

  if( !exists("eqtlInfo") || is.null(eqtlInfo)){
    message(i," | gene", traitGene, "has no eqtl associations, next!")
    message(" = None eQTL associations obtained of gene [",traitGene,"], please change the gene name or ENSEMBLE ID.")
    return(NULL)
  }
  # 如果有有多个 qtl_group, 只保留第一个。
  eqtlInfo <- eqtlInfo[ qtl_group ==eqtlInfo[1,]$qtl_group,]
  eqtlInfo[,position:=.(pos)]

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
  gwasDF <- gwasDF[,1:7]
  data.table::setDT(gwasDF)
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf", "beta", "se")
  # gwas subset:
  gwasDF <- na.omit(gwasDF)

  p_chrom_tmp <- ifelse(stringr::str_detect(gwasDF[1,]$chr, "chr"),P_chrom, stringr::str_remove(P_chrom, "chr"))
  gwasDF <- gwasDF[chr==p_chrom_tmp,]
  # convert variable class:
  gwasDF[,c("position", "pValue", "maf", "beta", "se"):=.(as.numeric(position), as.numeric(pValue), as.numeric(maf), as.numeric(beta), as.numeric(se))]
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
                                         gwasDF[,.(rsid, maf, pValue, beta, se)]
    )
    dataRanges_hg38 <- unlist(rtracklayer::liftOver(dataRanges, ch))
    # retain the SNPs located in chromosome 1:22
    dataRanges_hg38 <- dataRanges_hg38[which(as.character(GenomeInfoDb::seqnames(dataRanges_hg38)) %in% paste0("chr", 1:22)),]
    message("== ",length(dataRanges_hg38),"/",nrow(gwasDF)," (",round(length(dataRanges_hg38)/nrow(gwasDF)*100,2),"%)"," left.")
    gwasDFnew <- cbind(data.table::data.table(rsid = dataRanges_hg38$rsid, maf=dataRanges_hg38$maf, pValue=dataRanges_hg38$pValue, beta=dataRanges_hg38$beta, se=dataRanges_hg38$se),
                       data.table::data.table(chr = as.character(GenomeInfoDb::seqnames(dataRanges_hg38)), position=as.data.table(IRanges::ranges(dataRanges_hg38))$start ))
    gwasDF <- gwasDFnew[chr==P_chrom,.(rsid, chr, position, pValue, maf, beta, se)]
    rm(dataRanges, dataRanges_hg38, gwasDFnew)
  }


  #
  tissueSiteDetailId <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
  gwasEqtlInfo <- merge(gwasDF, eqtlInfo[,.(rsid, maf, pValue, position,beta, se)], by=c("rsid", "position"), suffixes = c(".gwas",".eqtl"))
  gwasEqtlInfo <- gwasEqtlInfo[se.eqtl !=0 & se.gwas !=0, ]

  if(nrow(gwasEqtlInfo)==0){
    message("No shared variants between eQTL and GWAS, please check your input!.")
    return(NULL)
  }
  message("== Start the colocalization analysis of gene [", traitGene,"]")

  if(method=="coloc"){
    message("== Using method: coloc")
    # 防止 check_dataset 中 p = pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2 出错
    suppressWarnings(coloc_Out <- coloc::coloc.abf(dataset1 = list( pvalues = gwasEqtlInfo$pValue.gwas, type="quant", N=gwasSampleNum, snp=gwasEqtlInfo$rsid, MAF=gwasEqtlInfo$maf.gwas),
                                                   dataset2 = list( pvalues = gwasEqtlInfo$pValue.eqtl, type="quant", N=ifelse(is.na(sampleNum[tissueSiteDetailId, on="tissueSiteDetailId"]$sampleNum), 300, sampleNum[tissueSiteDetailId, on="tissueSiteDetailId"]$sampleNum), snp=gwasEqtlInfo$rsid, MAF= gwasEqtlInfo$maf.eqtl)))
    coloc_Out_results <- as.data.table(coloc_Out$results)
    # coloc_Out_results$gene <- traitGenes[i]
    coloc_Out_summary <- as.data.table(t(as.data.frame(coloc_Out$summary)))
    coloc_Out_summary$traitGene <- traitGene
    coloc_Out_summary$candidate_snp <- coloc_Out_results[order(-SNP.PP.H4)][1,]$snp
    coloc_Out_summary$SNP.PP.H4 <- coloc_Out_results[order(-SNP.PP.H4)][1,]$SNP.PP.H4
    message("== Done")

    print(coloc_Out_summary)
    # coloc_Out_summary$pearsonCoor <- cor(-log(gwasEqtlInfo$pValue.gwas, 10),-log(gwasEqtlInfo$pValue.eqtl, 10), method = "pearson")

    return(list(coloc_Out_summary=coloc_Out_summary, gwasEqtlInfo=gwasEqtlInfo))
  }else if( method=="hyprcoloc" ){
    message("== Using method: hyprcoloc")
    # construct beta matrix and se matrix:
    betasMat <- as.matrix(gwasEqtlInfo[,.( eqtl=beta.eqtl, gwas=beta.gwas)])
    rownames(betasMat) <- gwasEqtlInfo$rsid
    sesMat <- as.matrix(gwasEqtlInfo[,.(eqtl=se.eqtl, gwas=se.gwas)])
    rownames(sesMat) <- gwasEqtlInfo$rsid
    res <- hyprcoloc::hyprcoloc(effect.est = betasMat, effect.se= sesMat, trait.names=colnames(betasMat), snp.id=rownames(betasMat), bb.alg=bb.alg);
    hyprcoloc_Out_summary <- as.data.table(res$results)
    hyprcoloc_Out_summary$traitGene <- traitGene
    hyprcoloc_Out_summary <- hyprcoloc_Out_summary[,.(traitGene, posterior_prob, regional_prob, candidate_snp, posterior_explained_by_snp)]
    message(hyprcoloc_Out_summary)
    message("== Done")

    return(list(hyprcoloc_Out_summary=hyprcoloc_Out_summary, gwasEqtlInfo=gwasEqtlInfo))
  }else if( method=="Both"){
    message("== Using methods: coloc AND hyprcoloc.")
    # coloc:
    suppressWarnings(coloc_Out <- coloc::coloc.abf(dataset1 = list( pvalues = gwasEqtlInfo$pValue.gwas, type="quant", N=gwasSampleNum, snp=gwasEqtlInfo$rsid, MAF=gwasEqtlInfo$maf.gwas),
                                                   dataset2 = list( pvalues = gwasEqtlInfo$pValue.eqtl, type="quant", N=ifelse(is.na(sampleNum[tissueSiteDetailId, on="tissueSiteDetailId"]$sampleNum), 300, sampleNum[tissueSiteDetailId, on="tissueSiteDetailId"]$sampleNum), snp=gwasEqtlInfo$rsid, MAF= gwasEqtlInfo$maf.eqtl)))
    coloc_Out_results <- as.data.table(coloc_Out$results)
    # coloc_Out_results$gene <- traitGenes[i]
    coloc_Out_summary <- as.data.table(t(as.data.frame(coloc_Out$summary)))
    coloc_Out_summary$traitGene <- traitGene
    coloc_Out_summary$candidate_snp <- coloc_Out_results[order(-SNP.PP.H4)][1,]$snp
    coloc_Out_summary$SNP.PP.H4 <- coloc_Out_results[order(-SNP.PP.H4)][1,]$SNP.PP.H4
    # print(coloc_Out_summary)

    message("")
    # hyprcoloc:
    betasMat <- as.matrix(gwasEqtlInfo[,.( eqtl=beta.eqtl, gwas=beta.gwas)])
    rownames(betasMat) <- gwasEqtlInfo$rsid
    sesMat <- as.matrix(gwasEqtlInfo[,.(eqtl=se.eqtl, gwas=se.gwas)])
    rownames(sesMat) <- gwasEqtlInfo$rsid
    res <- hyprcoloc::hyprcoloc(effect.est = betasMat, effect.se= sesMat, trait.names=colnames(betasMat), snp.id=rownames(betasMat), bb.alg=FALSE);
    hyprcoloc_Out_summary <- as.data.table(res$results)
    hyprcoloc_Out_summary$traitGene <- traitGene
    hyprcoloc_Out_summary <- hyprcoloc_Out_summary[,.(traitGene, hypr_posterior=posterior_prob, hypr_regional_prob=regional_prob, hypr_candidate_snp=candidate_snp, hypr_posterior_explainedBySnp=posterior_explained_by_snp)]
    colocOut <- merge(coloc_Out_summary, hyprcoloc_Out_summary, by="traitGene" )
    print(colocOut)

    message("== Done")
    return(list(colocOut=colocOut, gwasEqtlInfo=gwasEqtlInfo))
  }
}


#' @title conduct colocalization analysis with customized QTL data
#'
#' @param gwasDF data.frame or data.table, required cols: rsid, chrom, position, pValue, maf, beta, se
#' @param qtlDF data.frame or data.table, required cols:  rsid, chrom, position, pValue, maf, beta, se
#' @param mafThreshold 0.01
#' @param gwasSampleNum 50000
#' @param qtlSampleNum 10000
#' @param method "coloc", "hyprcoloc", and "Both"
#' @param bb.alg TRUE
#'
#' @return A list
#' @export
#' @examples
#' \donttest{
#' url1 <- "https://raw.githubusercontent.com/dingruofan/exampleData/master/gwasDFsub_MMP7.txt"
#' url2 <- "https://raw.githubusercontent.com/dingruofan/exampleData/master/eqtl/MMP7_qtlDF.txt"
#' gwasDF <- data.table::fread(url1)
#' qtlDF <- data.table::fread(url2)
#' output <- xQTLanalyze_coloc_local(gwasDF = gwasDF, qtlDF=qtlDF, method="hyprcoloc")
#' }
xQTLanalyze_coloc_diy <- function(gwasDF, qtlDF, mafThreshold=0.01, gwasSampleNum=50000, qtlSampleNum=10000, method="coloc", bb.alg=FALSE){
  rsid <- chrom <- chr <- position <- se <- pValue <- snpId <- maf <- pos <- i <- variantId <- se.eqtl <- se.gwas <-SNP.PP.H4 <- beta.eqtl <- beta.gwas <- posterior_prob <- regional_prob <-candidate_snp <- posterior_explained_by_snp  <- NULL
  . <- NULL

  if(method == "coloc" && !requireNamespace("coloc")){
    stop("please install package \"coloc\" with install.packages(\"coloc\").")
  }

  if(method == "hyprcoloc" && !requireNamespace("hyprcoloc")){
    stop("please install package \"hyprcoloc\".")
  }

  if(method == "Both"){
    if( !requireNamespace("coloc") ){ stop("please install package \"coloc\" with install.packages(\"coloc\").") }
    if( !requireNamespace("hyprcoloc") ){ stop("please install package \"hyprcoloc\".") }
  }


  # eqtl dataset:
  names(qtlDF) <- c("rsid", "chrom",  "position", "pValue", "maf", "beta", "se")
  # chromosome:
  P_chrom <- paste0("chr",qtlDF[1,]$chrom)

  eqtlInfo <- qtlDF[,.(rsid,chrom, maf=as.numeric(maf), beta=as.numeric(beta), se=as.numeric(se), pValue=as.numeric(pValue), position=as.numeric(position))]

  eqtlInfo<- eqtlInfo[maf>mafThreshold & maf <1,]
  # 去重：
  eqtlInfo <- eqtlInfo[order(rsid, pValue)][!duplicated(rsid)]

  if(nrow(eqtlInfo)==0){
    message("Number of eQTL association is <1")
    return(NULL)
  }

  message("Data processing...")

  ##################### gwas dataset:
  gwasDF <- gwasDF[,1:7]
  data.table::setDT(gwasDF)
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf", "beta", "se")
  # gwas subset:
  gwasDF <- na.omit(gwasDF)

  p_chrom_tmp <- ifelse(stringr::str_detect(gwasDF[1,]$chr, "chr"),P_chrom, stringr::str_remove(P_chrom, "chr"))
  gwasDF <- gwasDF[chr==p_chrom_tmp,]
  # convert variable class:
  gwasDF[,c("position", "pValue", "maf", "beta", "se"):=.(as.numeric(position), as.numeric(pValue), as.numeric(maf), as.numeric(beta), as.numeric(se))]
  # MAF filter:
  gwasDF <- gwasDF[maf > mafThreshold & maf<1,]
  # 去重：
  gwasDF <- gwasDF[order(rsid, pValue)][!duplicated(rsid)]
  # retain SNPs with rs id:
  gwasDF <- gwasDF[stringr::str_detect(rsid,stringr::regex("^rs")),]

  message(nrow(gwasDF))
  # chromosome revise：
  if( !stringr::str_detect(gwasDF[1,]$chr, stringr::regex("^chr")) ){
    gwasDF$chr <- paste0("chr", gwasDF$chr)
  }


  #
  gwasEqtlInfo <- merge(gwasDF, eqtlInfo[,.(rsid, maf, pValue, position,beta, se)], by=c("rsid", "position"), suffixes = c(".gwas",".eqtl"))
  gwasEqtlInfo <- gwasEqtlInfo[se.eqtl !=0 & se.gwas !=0, ]

  if(nrow(gwasEqtlInfo)==0){
    message("No shared variants between eQTL and GWAS, please check your input!.")
    return(NULL)
  }
  message("== Start the colocalization analysis")

  if(method=="coloc"){
    message("== Using method: coloc")
    # 防止 check_dataset 中 p = pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2 出错
    suppressWarnings(coloc_Out <- coloc::coloc.abf(dataset1 = list( pvalues = gwasEqtlInfo$pValue.gwas, type="quant", N=gwasSampleNum, snp=gwasEqtlInfo$rsid, MAF=gwasEqtlInfo$maf.gwas),
                                                   dataset2 = list( pvalues = gwasEqtlInfo$pValue.eqtl, type="quant", N=qtlSampleNum, snp=gwasEqtlInfo$rsid, MAF= gwasEqtlInfo$maf.eqtl)))
    coloc_Out_results <- as.data.table(coloc_Out$results)
    # coloc_Out_results$gene <- traitGenes[i]
    coloc_Out_summary <- as.data.table(t(as.data.frame(coloc_Out$summary)))
    coloc_Out_summary$candidate_snp <- coloc_Out_results[order(-SNP.PP.H4)][1,]$snp
    coloc_Out_summary$SNP.PP.H4 <- coloc_Out_results[order(-SNP.PP.H4)][1,]$SNP.PP.H4
    message("== Done")

    print(coloc_Out_summary)
    # coloc_Out_summary$pearsonCoor <- cor(-log(gwasEqtlInfo$pValue.gwas, 10),-log(gwasEqtlInfo$pValue.eqtl, 10), method = "pearson")

    return(list(coloc_Out_summary=coloc_Out_summary, gwasEqtlInfo=gwasEqtlInfo))
  }else if( method=="hyprcoloc" ){
    message("== Using method: hyprcoloc")
    # construct beta matrix and se matrix:
    betasMat <- as.matrix(gwasEqtlInfo[,.( eqtl=beta.eqtl, gwas=beta.gwas)])
    rownames(betasMat) <- gwasEqtlInfo$rsid
    sesMat <- as.matrix(gwasEqtlInfo[,.(eqtl=se.eqtl, gwas=se.gwas)])
    rownames(sesMat) <- gwasEqtlInfo$rsid
    res <- hyprcoloc::hyprcoloc(effect.est = betasMat, effect.se= sesMat, trait.names=colnames(betasMat), snp.id=rownames(betasMat), bb.alg=bb.alg);
    hyprcoloc_Out_summary <- as.data.table(res$results)
    hyprcoloc_Out_summary <- hyprcoloc_Out_summary[,.(posterior_prob, regional_prob, candidate_snp, posterior_explained_by_snp)]
    message(hyprcoloc_Out_summary)
    message("== Done")

    return(list(hyprcoloc_Out_summary=hyprcoloc_Out_summary, gwasEqtlInfo=gwasEqtlInfo))
  }else if( method=="Both"){
    message("== Using methods: coloc AND hyprcoloc.")
    # coloc:
    suppressWarnings(coloc_Out <- coloc::coloc.abf(dataset1 = list( pvalues = gwasEqtlInfo$pValue.gwas, type="quant", N=gwasSampleNum, snp=gwasEqtlInfo$rsid, MAF=gwasEqtlInfo$maf.gwas),
                                                   dataset2 = list( pvalues = gwasEqtlInfo$pValue.eqtl, type="quant", N=qtlSampleNum, snp=gwasEqtlInfo$rsid, MAF= gwasEqtlInfo$maf.eqtl)))
    coloc_Out_results <- as.data.table(coloc_Out$results)
    # coloc_Out_results$gene <- traitGenes[i]
    coloc_Out_summary <- as.data.table(t(as.data.frame(coloc_Out$summary)))
    coloc_Out_summary$candidate_snp <- coloc_Out_results[order(-SNP.PP.H4)][1,]$snp
    coloc_Out_summary$SNP.PP.H4 <- coloc_Out_results[order(-SNP.PP.H4)][1,]$SNP.PP.H4
    # print(coloc_Out_summary)

    message("")
    # hyprcoloc:
    betasMat <- as.matrix(gwasEqtlInfo[,.( eqtl=beta.eqtl, gwas=beta.gwas)])
    rownames(betasMat) <- gwasEqtlInfo$rsid
    sesMat <- as.matrix(gwasEqtlInfo[,.(eqtl=se.eqtl, gwas=se.gwas)])
    rownames(sesMat) <- gwasEqtlInfo$rsid
    res <- hyprcoloc::hyprcoloc(effect.est = betasMat, effect.se= sesMat, trait.names=colnames(betasMat), snp.id=rownames(betasMat), bb.alg=FALSE);
    hyprcoloc_Out_summary <- as.data.table(res$results)
    hyprcoloc_Out_summary <- hyprcoloc_Out_summary[,.( hypr_posterior=posterior_prob, hypr_regional_prob=regional_prob, hypr_candidate_snp=candidate_snp, hypr_posterior_explainedBySnp=posterior_explained_by_snp)]
    colocOut <- cbind(coloc_Out_summary, hyprcoloc_Out_summary )
    print(colocOut)

    message("== Done")
    return(list(colocOut=colocOut, gwasEqtlInfo=gwasEqtlInfo))
  }
}



#' @title eQTL-specific analysis
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueLabels (a character vector) can be listed with `ebi_study_tissues`. If is null, use all tissue / cell-types. (Default)
#' @param study (character) Studies can be listed using `ebi_study_tissues`. If is null, use all studies (Default).
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @import data.table
#' @import stringr
#' @return A list containing four data.table objects, including: "snpLD" for LD details of the specified SNP; "assoAllLd" for eQTL details of LD-associated SNPs;  "lm_R2_logP" for liner regression results; "cor_R2_logP" for correlation outputs;
#' @export
#'
#' @examples
#' \donttest{
#' propensityRes <- xQTLanalyze_propensity( gene="MMP7", variantName="rs11568818", study="TwinsUK")
#' # propensityRes <- xQTLanalyze_propensity(gene="FLOT1", variantName="rs3130356", study="")
#' xQTLvisual_qtlPropensity(propensityRes)
#' }
xQTLanalyze_propensity <- function(gene="", geneType="auto", variantName="", variantType="auto", tissueLabels="", study="", population="EUR"){
  .<- slope <- logP_minMax <- NULL
  study_accession <- tissue_label <- pos <- snpPanel <- variantId <- R2 <-snpId <- pValue <- tissue <- study_id <- qtl_group <- SNP_B <- LDbins <- logP <- corRP <- NULL
  pValue_propensity <- pValue_eQTL <- NULL

  # binNum A integer value. Number of bins to split values of R2 of LD into homogeneous bins.
  binNum=4
  ebi_ST <- data.table::copy(ebi_study_tissues)

  # study- tissue:
  if( all(tissueLabels!="") & study!=""){
    ebi_ST <- ebi_ST[which(tolower(tissue_label) %in% tolower(tissueLabels) & tolower(ebi_ST$study_accession) == tolower(study) ), ][,.SD[1,], by="tissue_label"]
    message(study)
  }

  # all tissue- study:
  if(study=="" & all(tissueLabels=="") ){
    ebi_ST <- ebi_ST[,.SD[1,], by="tissue_label"]
  }

  if(study!="" & all(tissueLabels=="") ){
    ebi_ST <- ebi_ST[ which(tolower(ebi_ST$study_accession) == tolower(study)),][,.SD[1,], by="tissue_label"]
  }

  if(study=="" & all(tissueLabels!="") ){
    ebi_ST <- ebi_ST[ which(tolower(tissue_label) %in% tolower(tissueLabels)),][,.SD[1,], by="tissue_label"]
  }

  if(nrow(ebi_ST)==0){
    stop("Please check study id or tissue label.")
  }

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
  if( (!exists("variantInfo")) || nrow(variantInfo)==0){
    stop("Variant [",variantName,"] is not exists in GTEx, please check your input.")
  }
  if(nrow(variantInfo)>1){
    variantInfo <- variantInfo[1,]
  }

  # fetch LD:
  message("== Retrieve LD information of SNP: [",variantInfo$snpId,"]...")
  try(snpLD <- retrieveLD(variantInfo$chromosome, variantInfo$snpId, population))
  data.table::setDT(snpLD)
  snpLD <- snpLD[which(stringr::str_detect(SNP_B, stringr::regex("^rs"))), ]
  # try(snpLD <- retrieveLD_LDproxy(targetSnp= variantInfo$snpId, population = population,  windowSize = 500000,genomeVersion = "grch38", token="9246d2db7917") )
  #
  if(nrow(snpLD)<1){
    stop("No LD found for variant: [", variantInfo$snpId, "]")
  }
  message("== Number of LD-associated variants: ", nrow(snpLD))

  # remove variants don't exist in GTEx.
  # chrom_snp <- stringr::str_split(snpLD[1,]$Coord, ":")[[1]][1]
  # snpLD$pos <- as.numeric(unlist(lapply(snpLD$Coord, function(x){ stringr::str_split(x,":")[[1]][2] })))
  # a<- xQTLquery_varPos(chrom = chrom_snp, pos = snpLD$pos, datasetId = "gtex_v8")

  # 从LD 0-1的10个区间随机选10个突变进行作图：
  # snpLD <- snpLD[,.(SNP_A= variantInfo$snpId, SNP_B=RS_Number, R2)]
  # # cut LD into 10 bins
  if( nrow(snpLD)>100){
    message("== Select by bins")
    binNumForSample_tmp <- 10
    varNumForSample_tmp <- 10
    snpLD$LDbins <- as.character(cut(snpLD$R2, breaks=seq(0,1,length.out=(binNumForSample_tmp+1)) ))
    snpLDForSample_tmp <- as.data.frame(tapply( 1:nrow(snpLD), snpLD$LDbins, function(x){ if(length(x)<=varNumForSample_tmp){return(x)}else{ return(sample(x, varNumForSample_tmp, replace = FALSE)) } }))
    names(snpLDForSample_tmp) <- "ID"
    snpLDForSample_tmp$LDbins <- rownames(snpLDForSample_tmp)
    data.table::setDT(snpLDForSample_tmp)
    # #
    snpLD <- snpLD[do.call(c,snpLDForSample_tmp$ID),-c("LDbins")][order(R2)]
    rm(binNumForSample_tmp, varNumForSample_tmp, snpLDForSample_tmp)
  }

  # Hypothesis testing
  # 对于每个组织分别进行假设检验
  # 1. 获取所有组织该基因的 eQTL. 2. 为每个 LD-associated gene构建panel. 3. 获得real panel和simulated panel 的SNP. 4.
  ebi_ST$pValue_eQTL <- (-1)
  ebi_ST$pValue_propensity <- (-1)
  assoAll <- data.table()
  for( i in 1:nrow(ebi_ST)){
    # 如果有多个 qtl group:
    message("")
    message("==> For tissue ",ebi_ST[i,]$tissue_label, " (",i,"/",nrow(ebi_ST),")")
    geneAsso_i <- xQTLdownload_eqtlAllAsso(gene=gene, geneType = geneType, tissueLabel = ebi_ST[i,]$tissue_label, study=ebi_ST[i,]$study_accession)
    if(!exists("geneAsso_i") ||is.null(geneAsso_i)|| nrow(geneAsso_i)==0){
      next()
    }else{
      geneAsso_i <- geneAsso_i[ qtl_group ==geneAsso_i[1,]$qtl_group,][order(pos)]
    }
    # LD-associated gene for plot:
    asso_I <- geneAsso_i[snpId %in% union(snpLD$SNP_B, snpLD$SNP_A),.(snpId, pValue, beta, tissue, tissue_label, study_id, qtl_group)]
    assoAll <- rbind(assoAll, asso_I)
    rm(asso_I)

    # all variants of gene:
    assoAllLd_i <- geneAsso_i[,.(snpId, pos, pValue)]
    assoAllLd_i <- merge(assoAllLd_i, snpLD[,.(snpId=SNP_B,R2)], by="snpId", sort=FALSE)
    assoAllLd_i$snpPanel <- paste0("s",1:nrow(assoAllLd_i))

    assoControl <- rbindlist(lapply(1:nrow(assoAllLd_i), function(snp_j){
      tmp <- geneAsso_i[data.table(ID=1:nrow(geneAsso_i), pos = geneAsso_i$pos, diff=abs(geneAsso_i$pos-assoAllLd_i[snp_j,]$pos))[diff!=0][order(diff)][1:5,]$ID,.(snpId, pos, pValue)]
      tmp$snpPanel <- assoAllLd_i[snp_j,]$snpPanel
      return(tmp)
    }))
    #

    asso_i <- rbind(assoAllLd_i[,-c("R2")], assoControl)
    asso_i <- merge(asso_i,  assoAllLd_i[,.(snpPanel,R2)], by="snpPanel")
    message("==> Start calculating p-value in ",ebi_ST[i,]$tissue_label, "; ",i,"/",nrow(ebi_ST))
    corValues <- unlist(lapply(1:1000, function(x){
      if(x %% 100 ==0){ message("==> Complete ",x/10,"% in tissue: ",  ebi_ST[i,]$tissue_label ," (",i,"/",nrow(ebi_ST),")")}
      sampleSnps <- rbindlist(lapply(unique(asso_i$snpPanel), function(xx){ a=asso_i[snpPanel==xx];a[sample(1:nrow(a),1),] }))
      return(cor(sampleSnps$R2, -log10(sampleSnps$pValue)))
    }))
    corReal <- cor(assoAllLd_i$R2, -log10(assoAllLd_i$pValue))
    ebi_ST[i, "pValue_propensity"] <- 1-(length(which(corReal>corValues)))/1001
    ebi_ST[i, "pValue_eQTL"] <- geneAsso_i[snpId==variantInfo$snpId,]$pValue
    rm(geneAsso_i, assoAllLd_i, assoControl, asso_i, corValues, corReal)
  }
  ebi_ST[pValue_propensity == (-1),"pValue_propensity"] <-NA
  ebi_ST[pValue_eQTL == (-1),"pValue_eQTL"] <-NA
  tissuePropensity  <- copy(ebi_ST)

  if(all(is.na(tissuePropensity$pValue_propensity))){
    stop("No associations were fetched from ", study)
  }

  # Retain min pvalue in each tissue_label-study, due to the duplication induced by qtl_group.
  assoAll <- assoAll[,.SD[which.min(pValue),], by=c("tissue", "tissue_label", "study_id", "qtl_group", "snpId")]
  assoAllLd <- merge(snpLD[,.(snpId=SNP_B, R2)], assoAll, by="snpId")
  assoAllLd$logP <- ifelse(assoAllLd$pValue==0, 0, (-log(assoAllLd$pValue,10)))
  # scale logP by tissue group:
  assoAllLd <- assoAllLd[,.(snpId, R2, beta, logP, logP_minMax=(logP-min(logP))/(max(logP)-min(logP))),by="tissue_label"]

  # saveRDS(ebi_ST,"../ebi_ST_FLOT1.rds")
  # query eQTL in each tissue:
  # eQTLs <- data.table()
  # for( i in 1:nrow(ebi_ST)){
  #   eQTL_i <- xQTLdownload_eqtlAllAsso(gene=gene, variantName = variantName, tissueLabel = ebi_ST[i,]$tissue_label, study=ebi_ST[i,]$study_id)
  #   if(!exists("eQTL_i") ||is.null(eQTL_i)|| nrow(eQTL_i)==0){
  #     next()
  #   }else{
  #     eQTL_i <- eQTL_i[,.(variantId, snpId, pos, tissue_label, study_id, pValue)]
  #   }
  #   eQTLs <- rbind(eQTLs, eQTL_i)
  # }
  # tissueSpecificity <- merge(ebi_ST, eQTLs[,.(tissue_label, study_id, pValue)], by=c("tissue_label", "study_id"))

  # lm:
  lm_f <- function(DT){
    lm_tmp <- lm(logP_minMax~R2,DT);
    return( data.table(slope=lm_tmp$coefficients['R2'], intercept= lm_tmp$coefficients['(Intercept)']) )}
  lm_R2_logP <- assoAllLd[,lm_f(.SD), by="tissue_label"][order(slope)]

  # cor:
  cor_R2_logP <- assoAllLd[,.(corRP=cor(R2, logP_minMax), corPvalue=cor.test(R2, logP_minMax)$p.value),by="tissue_label"][order(corRP)]
  cor_R2_logP$logCorP <- log10( cor_R2_logP$corPvalue)*(-1)

  return(list(snpLD=snpLD, tissuePropensity=tissuePropensity, cor_R2_logP=cor_R2_logP, lm_R2_logP=lm_R2_logP, assoAllLd=assoAllLd))
}



#' @title annotate variants with genome positon
#'
#' @param snpInfo A data.table/data.frame with two or three columns: chromosome and position.
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
#' snpHits <- xQTLanalyze_anno(snpInfo)
#' }
xQTLanalyze_anno <- function(snpInfo="", genomeVersion="hg38"){
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


#' @title variant enrichment analysis
#'
#' @param snpInfo A data.table/data.frame with two or three columns: chromosome and position.
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
#' enrichHits <- xQTLanalyze_enrich(snpInfo,enrichElement="Enhancer")
#' }
xQTLanalyze_enrich <- function(snpInfo="",  genomeVersion="hg38", enrichElement="Promoter", distLimit=1e6){
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
    message("      promoters...")
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
    message("      enhancer...")
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
    message("      tf cluster...")
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





#' @title calculate genomic control inflation factor
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
#' xQTLanalyze_calLambda(qtl[,.(pValue)])
#'
#' # calculate lambda value for each group:
#' qtl$groups <- sample(c(0,1),size = nrow(qtl), replace = TRUE)
#' xQTLanalyze_calLambda(qtl[,.(pValue, groups)])
#' }
xQTLanalyze_calLambda <- function(summaryDT){
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



