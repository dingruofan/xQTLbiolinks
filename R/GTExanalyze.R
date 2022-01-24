#' @title Sentinel SNPs detection in GWAS data.
#' @description Return sentinel snps whose pValue < 5e-8(default) and SNP-to-SNP distance > 1e6 bp.
#' @param gwasDF GWAS data.frame
#' @import data.table
#' @import stringr
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @importFrom  rtracklayer import.chain liftOver
#' @importFrom GenomeInfoDb seqnames
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#'    gwasFile <- tempfile(pattern = "file")
#'    gwasURL <- "https://raw.githubusercontent.com/dingruofan/exampleData/master/gwas/AD/GLGC_AD_chr1_6_Sub3.txt"
#'    utils::download.file(gwasURL, destfile=gwasFile)
#'    gwasDF <- data.table::fread(gwasFile, sep="\t")
#'    gwasDF <- gwasDF[, .(rsid, chr, position, P, maf)]
#'    sentinelSnpDF <- GTExanalyze_getSentinelSnp(gwasDF)
#'
#'    gwasDF <- fread("D:\\R_project\\GLGC_CG0052_result.txt.gz", sep="\t")
#'    gwasDF <- gwasDF[,.(rsid, chr, position, `p-value`, maf)]
#'    sentinelSnpDF_hg38 <- GTExanalyze_getSentinelSnp(gwasDF, centerRange=1e6, genomeVersion="grch37", grch37To38=TRUE)
#'    sentinelSnpDF_hg19 <- GTExanalyze_getSentinelSnp(gwasDF, centerRange=1e6, genomeVersion="grch37", grch37To38=FALSE)
#' }
GTExanalyze_getSentinelSnp <- function(gwasDF, pValueThreshold=5e-8, centerRange=1e4, mafThreshold = 0.01, genomeVersion="grch38", grch37To38 = FALSE){
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
    path = system.file(package="GTExbiolinks", "extdata", "hg19ToHg38.over.chain")
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

#' @title Detect the genes around the sentinel SNPs
#'
#' @param sentinelSnpDF A data.table. Better be the results from the function "GTExanalyze_getSentinelSnp", five columns are required, including "rsid", "chr", "position", "pValue", and "maf".
#' @param colocRange A integer value. Window size centered on SNP. Default: 1e6 bp
#' @param genomeVersion "grch38" or "grch19". Default: "grch38"
#' @import data.table
#' @import stringr
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @importFrom  rtracklayer import.chain liftOver
#' @importFrom GenomeInfoDb seqnames
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#'   sentinelSnpsURL <- "https://gitee.com/stronghoney/exampleData/raw/master/gwas/GLGC_CG0052/sentinelSnpDF.txt"
#'   sentinelSnpDF <- data.table::fread(rawToChar(curl::curl_fetch_memory(sentinelSnpsURL)$content), sep="\t")
#'   traitsAll_fromh19 <- GTExanalyze_getTraits(sentinelSnpDF_hg19, colocRange=1e6, genomeVersion="grch37", grch37To38=TRUE )
#'   traitsAll <- GTExanalyze_getTraits(sentinelSnpDF_hg38, colocRange=1e6, genomeVersion="grch38", grch37To38=FALSE )
#' }
GTExanalyze_getTraits <- function(sentinelSnpDF, colocRange=1e6, genomeVersion="grch38", grch37To38=FALSE){

  data.table::as.data.table(sentinelSnpDF)

  # (未做) 由于下一步的 GTExdownload_eqtlAllPost 函数只能 query 基于 hg38(v26) 的突变 1e6 bp附近的基因，所以如果输入的GWAS是 hg19 的突变坐标，需要进行转换为38，然后再进行下一步 eqtl sentinel snp filter.
  # 由于从 EBI category 里获得的是 hg38(v26) 的信息，所以如果这一步是 hg19 的1e6范围内，则在 hg38里就会未必，所以需要这一步，如果是hg19，则对突变的坐标进行变换：
  ####################### convert hg19 to hg38:
  if(genomeVersion =="grch37" & grch37To38){
    message("== Converting SNPs' coordinate to GRCH38... ")
    path = system.file(package="GTExbiolinks", "extdata", "hg19ToHg38.over.chain")
    ch = rtracklayer::import.chain(path)
    dataRanges <- GenomicRanges::GRanges(sentinelSnpDF$chr,
                                         IRanges::IRanges(sentinelSnpDF$position, sentinelSnpDF$position),
                                         strand= "*",
                                         sentinelSnpDF[,.(rsid, maf, pValue)]
    )
    dataRanges_hg38 <- unlist(rtracklayer::liftOver(dataRanges, ch))
    # retain the SNPs located in chromosome 1:22
    dataRanges_hg38 <- dataRanges_hg38[which(as.character(GenomeInfoDb::seqnames(dataRanges_hg38)) %in% paste0("chr", 1:22)),]
    sentinelSnpDF <- cbind(data.table(rsid = dataRanges_hg38$rsid, maf=dataRanges_hg38$maf, pValue=dataRanges_hg38$pValue),
                    data.table::data.table(chr = as.character(GenomeInfoDb::seqnames(dataRanges_hg38)), position=as.data.table(IRanges::ranges(dataRanges_hg38))$start ))
    sentinelSnpDF <- sentinelSnpDF[,.(rsid, chr, position, pValue, maf)]
    message("== ",length(dataRanges_hg38),"/",nrow(sentinelSnpDF)," left.")
    rm(dataRanges, dataRanges_hg38)
  }else if(genomeVersion =="grch38" & grch37To38){
    stop("Only grch37 genome version can be converted to grch38!")
  }else if(!genomeVersion %in% c("grch38", "grch37")){
    stop("Paramater genomeVersion must be choosen from grch38 and grch37, default: grch38.")
  }



  ###################### tissueSiteDetail="Brain - Cortex", maf_threshold=0.01
  geneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges, genomeVersion = "v26")
  datasetId="gtex_v8"
  GTExVersion="v26"

  chrAll <- unique(sentinelSnpDF$chr)
  chrAll <- chrAll[order(as.numeric(str_remove(chrAll,"chr")))]

  ########################## # 确保这些 SNPs 在 GTEx 里是有的， rsid 和 位置信息都一致才保留：
  # message("== Validating variant in GTEx....")
  # sentinelSnpDFNew <- data.table()
  # for(i in 1:length(chrAll)){
  #   sentinelSnpChromRaw <- sentinelSnpDF[chr==chrAll[i],]
  #   snpTMP <- suppressMessages(GTExquery_varPos(chrom = ifelse(datasetId=="gtex_v8", chrAll[i], str_remove(chrAll[i], "chr")), pos = sentinelSnpChromRaw$position, datasetId = datasetId))
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
    Traits <- rbindlist(lapply(1:nrow(sentinelSnpChrom), function(x){
      tmp<-sentinelSnpChrom[x,];
      startPos <- tmp$position - colocRange; endPos <- tmp$position + colocRange;
      traits <- geneInfoChrom[ (startPos>start & startPos<=end) | (endPos>start & endPos<=end) ];
      # traits$rsid <- tmp$rsid; traits$position <- tmp$position;traits$position <- tmp$position;
      traits[,c("rsid","position","pValue","maf"):=.(tmp$rsid, tmp$position, tmp$pValue, tmp$maf)];
      return(traits)
    }))
    Traits$chromosome <- chrAll[i]
    traitsAll <- rbind(traitsAll, Traits)
    message(i," - ", chrAll[i]," - ",nrow(Traits))
    rm(sentinelSnpChrom, geneInfoChrom, Traits)
  }
  message("== Fetching [", nrow(traitsAll),"] genes' information from the GTEx.")
  geneAnnot <- GTExquery_gene(unique(traitsAll$gencodeId), geneType = "gencodeId", gencodeVersion = GTExVersion)
  if( exists("geneAnnot") && !is.null(geneAnnot) && nrow(geneAnnot)>0 ){
    traitsAll <- merge(geneAnnot[,.(genes, geneSymbol,gencodeId, geneType, description)],traitsAll, by.x="genes",by.y="gencodeId",  all.x=TRUE)[,-c("genes")]
    traitsAll <- traitsAll[,.( chromosome, geneStart=start, geneEnd=end, geneStrand=strand, geneSymbol, gencodeId, rsid, position, pValue, maf )][order(as.numeric(str_remove(chromosome, "chr")), pValue, position)]
    message("== Totally, [",nrow(traitsAll), "] associations between [",length(unique(traitsAll$gencodeId)),"] traits genes and [",length(unique(traitsAll$rsid)),"] SNPs are detected." )
    return(traitsAll)
  }else{
    return(traitsAll)
  }

}

# check existence of association of sentinel snp - gene of taritAll.
# for(i in 1:nrow(a)){
#   # GTExdownload_eqtlAll(variantName = traitsAll[i,]$rsid, gene = traitsAll[i,]$gencodeId, geneType = "gencodeId", tissueSiteDetail = tissueSiteDetail)
#   aa <- suppressMessages( GTExdownload_eqtlAll(variantName =a[i]$snpId, gene = a[i]$gencodeId, geneType = "gencodeId", tissueSiteDetail = tissueSiteDetail) )
#   if( nrow(aa) >0){
#     message(i," 有数据")
#   }else{
#     message(i, " 无数据")
#   }
# }


#' @title Colocalization analysis with deteched trait
#'
#' @param gwasDF 1
#' @param traitGenes 1
#' @param tissueSiteDetail 1
#' @param population 1
#' @param token 1
#' @param method 1
#'
#' @return coloc resut
#' @export
#'
#' @examples
#' \donttest{
#'
#'   traitsAll <- merge(traitsAllEqtl[pValue !="N/A",.(gencodeId, rsid=snpId)], traitsAll, by = c("gencodeId", "rsid"))
#'   traitsAll <- traitsAll[order(pValue)]
#'
#'   traitsAllURL <- "https://gitee.com/stronghoney/exampleData/raw/master/gwas/AD/traitsAll.txt"
#'   traitsAll <- data.table::fread(rawToChar(curl::curl_fetch_memory(traitsAllURL)$content), sep="\t")

#'
#'   colocResultAll <- list()
#'   for(i in 154:nrow(traitsAll)){
#'     message(i," - ", date())
#'     traitGene <- traitsAll[order(pValue)][i,]$geneSymbol
#'     sentinelSnp <- traitsAll[order(pValue)][i,]$rsid
#'     colocResult <- GTExanalyze_coloc(gwasDF, traitGene=traitGene, sentinelSnp= sentinelSnp, tissueSiteDetail="Breast - Mammary Tissue")
#'     colocResultAll[[i]] <- colocResult
#'   }
#'   range( unlist(lapply(1:length(colocResultAll), function(x){colocResultAll[[x]]$coloc_Out_summary$PP.H4.abf})) )
#'
#'   17
#'   # 为啥全为 NA：
#'   a <- GTExdownload_eqtlAll(gene="ENSG00000152670.18", tissueSiteDetail = tissueSiteDetail)
#'   a1 <- GTExdownload_eqtlSig(gene="ENSG00000152670.18", tissueSiteDetail = tissueSiteDetail)
#'
#'   colocResult$gwasEqtlInfo
#'   # eQTL locuszoom:
#'   GTExvisual_locusZoom( colocResult$gwasEqtlInfo[,c("rsid","chr","position","pValue.eqtl")], population="EUR",genomeVersion="grch38" )
#'   # GWAS locuszoom:
#'   GTExvisual_locusZoom( colocResult$gwasEqtlInfo[,c("rsid","chr","position","pValue.gwas")], population="EUR",genomeVersion="grch38" )
#'   # locuscompare:
#'   GTExvisual_locusCompare( colocResult$gwasEqtlInfo[,c("rsid","pValue.eqtl")], colocResult$gwasEqtlInfo[,c("rsid","pValue.gwas")] )
#' }
GTExanalyze_coloc <- function(gwasDF, traitGene, geneType="geneSymbol", sentinelSnp, variantType="snpId", tissueSiteDetail="", colocRange=1e6, mafThreshold=0.01, gwasSampleNum=50000, population="EUR", token= "9246d2db7917", method="coloc"){
  # tissueSiteDetail="Brain - Cortex"
  # colocRange=1e6
  # numSnp = 300
  # gwasSampleNum=50000


  ###################### eqtl dataset:
  eqtlInfo <- GTExdownload_assoAll(traitGene,geneType = geneType, tissueSiteDetail=tissueSiteDetail, withdbSNPID = FALSE)
  if( is.null(eqtlInfo) || nrow(eqtlInfo)==0){
    message(i," | gene", traitGene, "has no eqtl associations, next!")
    message(" = None eQTL associations obtained of gene [",traitGene,"], please change the gene name or ENSEMBLE ID.")
    return(NULL)
  }
  eqtlInfo[,"position":= .( lapply(variantId, function(x){str_split(x, fixed("_"))[[1]][2]}) )]
  eqtlInfo<- eqtlInfo[maf>mafThreshold & maf <1,]
  # 去重：
  eqtlInfo <- eqtlInfo[order(snpId, pValue)][!duplicated(snpId)]

  # chromosome:
  P_chrom <- str_split(eqtlInfo[1,]$variantId, fixed("_"))[[1]][1]

  eqtlInfo <- na.omit( eqtlInfo[,.(rsid=snpId, maf, beta, se, pValue)] )




  ##################### gwas dataset:
  gwasDF <- gwasDF[,1:5]
  data.table::setDT(gwasDF)
  message("== Start the colocalization analysis of gene ", traitGene)
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf")
  # gwas subset:
  gwasDF <- na.omit(gwasDF)
  # chromosome revise：
  if( !str_detect(gwasDF[1,]$chr, regex("^chr")) ){
    gwasDF$chr <- paste0("chr", gwasDF$chr)
  }
  gwasDF <- gwasDF[chr==P_chrom,]
  # convert variable class:
  gwasDF[,c("position", "pValue", "maf")] <- gwasDF[,.(position=as.numeric(position), pValue=as.numeric(pValue), maf=as.numeric(maf))]
  # MAF filter:
  gwasDF <- gwasDF[maf > mafThreshold & maf<1,]
  # 去重：
  gwasDF <- gwasDF[order(rsid, pValue)][!duplicated(rsid)]
  # retain SNPs with rs id:
  gwasDF <- gwasDF[stringr::str_detect(rsid,stringr::regex("^rs")),]


  #####################  2
  tissueSiteDetailId <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
  gwasEqtldata <- merge(gwasDF, eqtlInfo[,.(rsid, maf, pValue)], by="rsid", suffixes = c(".gwas",".eqtl"))
  # centerSnp <- gwasEqtldata[which.min(gwasEqtldata$pValue.gwas),]
  centerSnp <- GTExquery_varId(sentinelSnp, variantType = variantType)
  if(nrow(centerSnp)==0 ){
    centerSnp <- gwasDF[which.min(gwasDF$pValue),]
    gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$pos-(colocRange)) & position<=(centerSnp$pos+(colocRange)),][order(position)]
  }
  if( colocRange==0 ){
    gwasEqtlInfo <- gwasEqtldata[order(position)]
  }else {
    gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$pos-(colocRange)) & position<=(centerSnp$pos+(colocRange)),][order(position)]
  }
  if(nrow(gwasEqtlInfo)==0){
    centerSnp <- gwasDF[which.min(gwasDF$pValue),]
    gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$pos-(colocRange)) & position<=(centerSnp$pos+(colocRange)),][order(position)]
  }
  if(nrow(gwasEqtlInfo)==0){
    return(NULL)
  }

  # 只选择中心附近的 Num snp 进行分析：
  # 防止 check_dataset中 p = pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2 出错
  suppressWarnings(coloc_Out <- coloc::coloc.abf(dataset1 = list( pvalues = gwasEqtlInfo$pValue.gwas, type="quant", N=gwasSampleNum, snp=gwasEqtlInfo$rsid, MAF=gwasEqtlInfo$maf.gwas),
                                                 dataset2 = list( pvalues = gwasEqtlInfo$pValue.eqtl, type="quant", N=sampleNum[tissueSiteDetailId, on="tissueSiteDetailId"]$sampleNum, snp=gwasEqtlInfo$rsid, MAF= gwasEqtlInfo$maf.eqtl)))
  coloc_Out_results <- as.data.table(coloc_Out$results)
  # coloc_Out_results$gene <- traitGenes[i]
  coloc_Out_summary <- as.data.table(t(as.data.frame(coloc_Out$summary)))
  # coloc_Out_summary$pearsonCoor <- cor(-log(gwasEqtlInfo$pValue.gwas, 10),-log(gwasEqtlInfo$pValue.eqtl, 10), method = "pearson")

  return(list(coloc_Out_summary=coloc_Out_summary,coloc_Out_results=coloc_Out_results, gwasEqtlInfo=gwasEqtlInfo))
}
