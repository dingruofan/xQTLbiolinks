#' @title Sentinel SNPs detection in GWAS data.
#' @description Return sentinel snps whose pValue < 5e-8(default) and SNP-to-SNP distance > 1e6 bp.
#' @param gwasDF GWAS data.frame
#' @import data.table
#' @import stringr
#' @import LDlinkR
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import GenomicFeatures
#' @import org.Hs.eg.db
#' @import GenomeInfoDb
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
#'    sentinelSnpDF <- GTExanalyze_getSentinelSnp(gwasDF, centerRange=1e4)
#' }
GTExanalyze_getSentinelSnp <- function(gwasDF, pValueThreshold=5e-8, centerRange=1e6, mafThreshold = 0.01){
  # Detect gene with sentinal SNP:
  gwasDF <- gwasDF[,1:5]
  message("== Start the detection of sentinel SNPs: ")
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

  # sentinel snp:
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
      startPos <- tmp[1,]$position-centerRange/2
      endPos <- tmp[1,]$position+centerRange/2
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
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#'   sentinelSnpsURL <- "https://gitee.com/stronghoney/exampleData/raw/master/gwas/GLGC_CG0052/sentinelSnpDF.txt"
#'   sentinelSnpDF <- data.table::fread(rawToChar(curl::curl_fetch_memory(sentinelSnpsURL)$content), sep="\t")
#'   traitsAll <- GTExanalyze_getTraits(sentinelSnpDF, colocRange=1e6, genomeVersion="grch37" )
#' }
GTExanalyze_getTraits <- function(sentinelSnpDF, colocRange=1e6, genomeVersion="grch38"){

  data.table::as.data.table(sentinelSnpDF)
  # tissueSiteDetail="Brain - Cortex", maf_threshold=0.01
  if(tolower(genomeVersion) == "grch38"){
    geneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges, genomeVersion = "v26")
    datasetId="gtex_v8"
  }else{
    geneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges, genomeVersion = "v19")
    datasetId = "gtex_v7"
  }


  #########################
  # 确保这些 SNPs 在 GTEx 里是有的， rsid 和 位置信息都一致才保留：
  message("== start variant validation.")
  chrAll <- unique(sentinelSnpDF$chr)
  chrAll <- chrAll[order(as.numeric(str_remove(chrAll,"chr")))]
  sentinelSnpDFNew <- data.table()
  for(i in 1:length(chrAll)){
    sentinelSnpChromRaw <- sentinelSnpDF[chr==chrAll[i],]
    snpTMP <- suppressMessages(GTExquery_varPos(chrom = ifelse(datasetId=="gtex_v8", chrAll[i], str_remove(chrAll[i], "chr")), pos = sentinelSnpChromRaw$position, datasetId = datasetId))
    # 使用 位置进行merge:， 因为 rsid 有可能有版本差异：
    sentinelSnpChrom <- merge(sentinelSnpChromRaw, unique(snpTMP[,.(position=pos,rsid = snpId)]), by=c("rsid", "position"))
    if(nrow(sentinelSnpChrom)>0){
      sentinelSnpDFNew <- rbind(sentinelSnpDFNew, sentinelSnpChrom)
    }else{
      next()
    }
    message("== ", chrAll[i], ", ", nrow(sentinelSnpChrom),"/", nrow(sentinelSnpChromRaw) ," retained in GTEx.")
    rm(sentinelSnpChromRaw, sentinelSnpChrom, snpTMP)
  }
  if(nrow(sentinelSnpDFNew)==0){
    stop("There is no shared variants found in your gwas dataset.")
  }else{
    message("== Totally, ",nrow(sentinelSnpDFNew), "/",nrow(sentinelSnpDF), " sentinel SNPs retained.")
  }

  ########################
  message("== Start the trait gene detection of sentinel snps. ")
  traitsAll <- data.table()
  for(i in 1:length(chrAll)){
    geneInfoChrom <- geneInfo[chromosome == chrAll[i]]
    sentinelSnpChrom <- sentinelSnpDFNew[chr==chrAll[i],]
    Traits <- rbindlist(lapply(1:nrow(sentinelSnpChrom), function(x){
      tmp<-sentinelSnpChrom[x,];
      startPos <- tmp$position - colocRange/2; endPos <- tmp$position + colocRange/2;
      traits <- geneInfoChrom[ (startPos>=start & startPos<=end) | (endPos>=start & endPos<=end) ];
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
  gencodeVersion <- ifelse(tolower(genomeVersion) == "grch38", "v26", "v19")
  geneAnnot <- GTExquery_gene(unique(traitsAll$gencodeId), geneType = "gencodeId", gencodeVersion = gencodeVersion)
  traitsAll <- merge(geneAnnot[,.(genes, geneSymbol,gencodeId, geneType, description)],traitsAll, by.x="genes",by.y="gencodeId",  all.x=TRUE)[,-c("genes")]

  message("== Totally, [",nrow(traitsAll), "] associations between [",length(unique(traitsAll$gencodeId)),"] traits genes and [",length(unique(traitsAll$rsid)),"] SNPs are detected." )
  return(traitsAll)
}

#' @title
#'
#' @param traitGenes
#' @param sentinelSnps
#' @param tissueSiteDetail
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'   traitGenes <- traitsAll$gencodeId
#'   sentinelSnps <- traitsAll$rsid
#' }
GTExanalyze_traitEqtlSig <- function(traitGenes, sentinelSnps, tissueSiteDetail=""){
  if(length(traitGenes)!= length(sentinelSnps)){
    stop("Number of traitgenes is not equal to that of sentinel SNPs!")
  }
  gtexEqtl <- GTExdownload_eqtlAllPost( genes= traitGenes, variants= sentinelSnps, tissueSiteDetail="Skin - Sun Exposed (Lower leg)", recordPerChunk = 100)
  return(gtexEqtl)
}


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
#'   traitsAllURL <- "https://gitee.com/stronghoney/exampleData/raw/master/gwas/AD/traitsAll.txt"
#'   traitsAll <- data.table::fread(rawToChar(curl::curl_fetch_memory(traitsAllURL)$content), sep="\t")
#'   traitGene <- traitsAll[3,]$geneSymbol
#'   sentinelSnp <- traitsAll[3,]$rsid
#'
#'   colocResult <- GTExanalyze_coloc(gwasDF, traitGenes)
#'   colocResult$gwasEqtlInfo
#'   # eQTL locuszoom:
#'   GTExvisual_locusZoom( colocResult$gwasEqtlInfo[,c("rsid","chr","position","pValue.eqtl")], population="EUR",genomeVersion="grch38" )
#'   # GWAS locuszoom:
#'   GTExvisual_locusZoom( colocResult$gwasEqtlInfo[,c("rsid","chr","position","pValue.gwas")], population="EUR",genomeVersion="grch38" )
#'   # locuscompare:
#'   GTExvisual_locusCompare( colocResult$gwasEqtlInfo[,c("rsid","pValue.eqtl")], colocResult$gwasEqtlInfo[,c("rsid","pValue.gwas")] )
#' }
GTExanalyze_coloc <- function(gwasDF, traitGene, sentinelSnp, tissueSiteDetail="Brain - Cortex", colocRange=1e6, mafThreshold=0.01, gwasSampleNum=50000, population="EUR", token= "9246d2db7917", method="coloc"){
  # tissueSiteDetail="Brain - Cortex"
  # colocRange=1e6
  # numSnp = 300
  # gwasSampleNum=50000

  #####################
  # eqtl dataset:
  eqtlInfo <- GTExdownload_assoAll(traitGene, tissueSiteDetail=tissueSiteDetail, withdbSNPID = FALSE)
  if(nrow(eqtlInfo)==0){
    message(i," | gene", traitGene, "has no eqtl associations, next!")
    stop(" = None eQTL associations obtained of gene [",traitGene,"], please change the gene name or ENSEMBLE ID.")
  }
  eqtlInfo[,"position":= .( lapply(variantId, function(x){str_split(x, fixed("_"))[[1]][2]}) )]
  eqtlInfo<- eqtlInfo[maf>mafThreshold & maf <1,]
  # 去重：
  eqtlInfo <- eqtlInfo[order(snpId, pValue)][!duplicated(snpId)]

  # chromosome:
  P_chrom <- str_split(eqtlInfo[1,]$variantId, fixed("_"))[[1]][1]

  eqtlInfo <- na.omit( eqtlInfo[,.(rsid=snpId, maf, beta, se, pValue)] )
  #####################


  #####################
  # gwas dataset:
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
  #####################

  #
  tissueSiteDetailId <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
  gwasEqtldata <- merge(gwasDF, eqtlInfo[,.(rsid, maf, pValue)], by="rsid", suffixes = c(".gwas",".eqtl"))
  # centerSnp <- gwasEqtldata[which.min(gwasEqtldata$pValue.gwas),]
  centerSnp <- gwasEqtldata[rsid == sentinelSnp,]
  if( nrow(centerSnp)==0){
    stop("None shared SNPs detected in eQTL and GWAS dataset, please choose new gene or SNPs.")
  }
  if( colocRange==0 ){
    gwasEqtlInfo <- gwasEqtldata[order(position)]
  }else {
    gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$position-(colocRange)/2) & position<=(centerSnp$position+(colocRange)/2),][order(position)]
  }
  # 只选择中心附近的 Num snp 进行分析：
  # minGwasPvalueVar <- which.min(gwasEqtlInfo$pValue.gwas)
  # gwasEqtlInfo <- gwasEqtlInfo[ifelse(minGwasPvalueVar>numSnp,minGwasPvalueVar-numSnp,1):ifelse( (nrow(gwasEqtlInfo)-minGwasPvalueVar)>numSnp, minGwasPvalueVar+numSnp, nrow(gwasEqtlInfo) ),]
  # 防止 check_dataset中 p = pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2 出错
  suppressWarnings(coloc_Out <- coloc::coloc.abf(dataset1 = list( pvalues = gwasEqtlInfo$pValue.gwas, type="quant", N=gwasSampleNum, snp=gwasEqtlInfo$rsid, MAF=gwasEqtlInfo$maf.gwas),
                                                 dataset2 = list( pvalues = gwasEqtlInfo$pValue.eqtl, type="quant", N=sampleNum[tissueSiteDetailId, on="tissueSiteDetailId"]$sampleNum, snp=gwasEqtlInfo$rsid, MAF= gwasEqtlInfo$maf.eqtl)))
  # coloc_Out_results <- as.data.table(coloc_Out$results)
  # coloc_Out_results$gene <- traitGenes[i]
  coloc_Out_summary <- as.data.table(t(as.data.frame(coloc_Out$summary)))
  coloc_Out_summary$pearsonCoor <- cor(-log(gwasEqtlInfo$pValue.gwas, 10),-log(gwasEqtlInfo$pValue.eqtl, 10), method = "pearson")
  coloc_Out_summary$gene <- traitGenes[i]
  if(exists("coloc_Out_summary") && nrow(coloc_Out_summary)>0){
    coloc_Out_summaryAll <- rbind(coloc_Out_summaryAll, coloc_Out_summary)
  }

  if( !exists("eqtlInfo") ){
    message("Note: No eQTL associations found for gene: ", paste0(traitGenes, collapse = ", "))
    return(NULL)
  }else{
    return(list(coloc_Out_summaryAll=coloc_Out_summaryAll, gwasEqtlInfo=gwasEqtlInfo))
  }
}
