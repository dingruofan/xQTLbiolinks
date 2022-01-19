#' @title Sentinel SNPs detection in GWAS data.
#'
#' @param gwasDF GWAS data.frame
#' @import data.table
#' @import stringr
#' @import LDlinkR
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import GenomicFeatures
#' @import org.Hs.eg.db
#' @import GenomeInfoDb
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'    gwasFile <- tempfile(pattern = "file")
#'    gwasURL <- "https://raw.githubusercontent.com/dingruofan/exampleData/master/gwas/AD/GLGC_AD_chr1_6_Sub3.txt"
#'    utils::download.file(gwasURL, destfile=gwasFile)
#'    gwasDF <- data.table::fread(gwasFile, sep="\t")
#'    gwasDF <- gwasDF[, .(rsid, chr, position, P, maf)]
#'    sentinelSnps <- GTExanalyze_getSentinelSnp(gwasDF)
#'
#'    gwasDF <- fread("D:\\R_project\\GLGC_CG0052_result.txt.gz", sep="\t")
#'    gwasDF <- gwasDF[,.(rsid, chr, position, `p-value`, maf)]
#'    sentinelSnps <- GTExanalyze_getSentinelSnp(gwasDF)
#' }
GTExanalyze_getSentinelSnp <- function(gwasDF, pValueThreshold=5e-8, centerRange=1e6 ){
  # Detect gene with sentinal SNP:
  gwasDF <- gwasDF[,1:5]
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf")
  gwasDF <- na.omit(gwasDF)
  # convert variable class:
  gwasDF[,c("position", "pValue", "maf")] <- gwasDF[,.(position=as.numeric(position), pValue=as.numeric(pValue), maf=as.numeric(maf))]

  # MAF filter:
  # gwasDF <- gwasDF[maf>maf_threshold & maf<1,]
  # 去重：
  # gwasDF <- gwasDF[order(rsid, pValue)][!duplicated(rsid)]
  # retain SNPs with rs id:
  # gwasDF <- gwasDF[stringr::str_detect(rsid,stringr::regex("^rs")),]

  # sentinel snp:
  gwasDFsub <- gwasDF[pValue<pValueThreshold, ][order(pValue, position)]

  chrAll <- unique(gwasDFsub$chr)
  sentinelSnps <- data.table()
  message("== Start the detection of sentinel SNPs: ")
  for(i in 1:length(chrAll)){
    gwasDFsubChrom <- gwasDFsub[chr==chrAll[i],]
    tmp <- copy(gwasDFsubChrom)
    sentinelSnps_count <- 0
    while(nrow(tmp)>0){
      sentinelSnps_count <- sentinelSnps_count+1
      sentinelSnps <- rbind(sentinelSnps, tmp[1,])
      startPos <- tmp[1,]$position-centerRange/2
      endPos <- tmp[1,]$position+centerRange/2
      tmp <- tmp[position<startPos | position>endPos][order(pValue, position)]
    }
    message("   In ",chrAll[i], ", ",sentinelSnps_count, " sentinel SNPs detected. ")
    rm(tmp, sentinelSnps_count, startPos, endPos)
  }
  message("== Totally, ", nrow(sentinelSnps), " sentinel SNPs located in ",length(unique(sentinelSnps$chr))  ," chromosomes have been detected!")
  return(sentinelSnps)
}

GTExanalyze_getTraits <- function(sentinelSnpDF, colocRange=1e6, genome="grch38"){
  # tissueSiteDetail="Brain - Cortex", maf_threshold=0.01
  if(genome == "grch38"){
    geneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges, genomeVersion = "V26")
  }else{
    geneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges, genomeVersion = "V19")
  }

  message("== Start : ")
  chrAll <- unique(sentinelSnpDF$chr)
  traitsAll <- data.table()
  for(i in 1:length(chrAll)){
    sentinelSnpChrom <- sentinelSnpDF[chr==chrAll[i],]
    geneInfoChrom <- geneInfo[chromosome == chrAll[i]]

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
  message("== End : ")
}


# 1. GTEx significant eQTL:
# gwasSNPeqtl <- rbindlist(lapply(1:nrow(gwasDFsub), function(i){
#   message(i,"/",nrow(gwasDFsub)," | " ,gwasDFsub$rsid[i])
#   tmp1 <- suppressMessages( GTExdownload_eqtlSig(variantName= gwasDFsub$rsid[i],  datasetId="gtex_v8", tissueSiteDetail="Brain - Cortex") )
#   if( !is.null(tmp1) && nrow(tmp1)!=0 ){
#     return(tmp1)
#   }else{
#     return(NULL)
#   }
# }))

# #2. Gene range:
# genes_txdb <- suppressMessages( GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene) )
# genes_txdb <- cbind(data.table(chr=as.character(GenomeInfoDb::seqnames(genes_txdb))), as.data.table(IRanges::ranges(genes_txdb)))
# genes_txdb <- genes_txdb[chr %in% paste0("chr",c(1:22,"X","Y"))]
# genes_txdb$start <-genes_txdb$start - upDownStream
# genes_txdb$end <-genes_txdb$end + upDownStream
# # merge:
# targetGene <- data.table()
# for( c in unique(gwasDFsub$chr)){
#   geneTmp <- genes_txdb[chr==c]
#   gwasDFsubTmp <- gwasDFsub[chr==c]
#   a <- data.table::rbindlist(lapply(1:nrow(gwasDFsubTmp), function(x){ pos1=as.integer(gwasDFsubTmp[x]$position); hited <- geneTmp[pos1 - start >=0 & pos1-end<=0 ]; if(nrow(hited)==0){return(NULL)}else{ hited <- cbind(hited, gwasDFsubTmp[x,c("rsid", "position")]);hited$id=x;return(hited) } }))
#   targetGene <- rbind(targetGene, a)
#   message( "chrom: ", c)
# }
# if(nrow(targetGene)==0){
#   message("Please reset pValueThreshold a greater value.")
#   return(NULL)
# }
# targetGeneUniq <- targetGene[,.( rsids=paste0(rsid, collapse = ", ")),by = c("chr", "start", "end", "width", "names")]
# # gene annotation:
# targetGeneAnno <- AnnotationDbi::select(org.Hs.eg.db, keys=targetGeneUniq$names, columns=c("SYMBOL", "ENSEMBL"), keytype="ENTREZID")
# targetGeneUniqAnno <- merge( targetGeneAnno, targetGeneUniq, by.x="ENTREZID", by.y="names")
# return(targetGeneUniqAnno)

# 3. ldlink
# tissueSiteDetailId <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
# a <- LDexpress(gwasDFsub$rsid[1:10],r2d_threshold=1, tissue = tissueSiteDetailId, pop = population, token = "9246d2db7917" )

#' @title coloc analysis with deteched trait
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
#'   traitGenes <- c("CR1")
#'   colocResult <- GTExanalyze_coloc(gwasDF, traitGenes)
#'   colocResult$gwasEqtlInfo
#'   # eQTL locuszoom:
#'   GTExvisual_locusZoom( colocResult$gwasEqtlInfo[,c("rsid","chr","position","pValue.eqtl")], population="EUR",genome="grch38" )
#'   # GWAS locuszoom:
#'   GTExvisual_locusZoom( colocResult$gwasEqtlInfo[,c("rsid","chr","position","pValue.gwas")], population="EUR",genome="grch38" )
#'   # locuscompare:
#'   GTExvisual_locusCompare( colocResult$gwasEqtlInfo[,c("rsid","pValue.eqtl")], colocResult$gwasEqtlInfo[,c("rsid","pValue.gwas")] )
#' }
GTExanalyze_coloc <- function(gwasDF, traitGenes, tissueSiteDetail="Brain - Cortex", colocRange=1e6, numSnp = 3000, gwasSampleNum=50000, population="EUR", token= "9246d2db7917", method="coloc"){
  # tissueSiteDetail="Brain - Cortex"
  # colocRange=1e6
  # numSnp = 300
  # gwasSampleNum=50000
  #
  gwasDF <- gwasDF[,1:5]
  names(gwasDF) <- c("rsid", "chr", "position", "pValue", "maf")
  gwasDF <- na.omit(gwasDF)
  # 去重：
  gwasDF <- gwasDF[order(rsid, pValue)][!duplicated(rsid)]
  #
  tissueSiteDetailId <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
  coloc_Out_summaryAll <- data.table()
  for( i in 1:length(traitGenes)){
    eqtlInfo <- GTExdownload_assoAll(traitGenes[i], tissueSiteDetail=tissueSiteDetail)
    if(nrow(eqtlInfo)==0){
      message(i," | gene", traitGenes[i], "has no eqtl associations, next!")
      next()
    }
    eqtlInfo<- eqtlInfo[maf>0,]
    # 去重：
    eqtlInfo <- eqtlInfo[!duplicated(snpId)]
    gwasEqtldata <- merge(gwasDF, eqtlInfo[,.(rsid=snpId, maf, beta, se, pValue)], by="rsid", suffixes = c(".gwas",".eqtl"))
    centerSnp <- gwasEqtldata[which.min(gwasEqtldata$pValue.gwas),]
    if( colocRange==0 ){
      gwasEqtlInfo <- gwasEqtldata[order(position)]
    }else{
      gwasEqtlInfo <- gwasEqtldata[position>=(centerSnp$position-(colocRange)/2) & position<=(centerSnp$position+(colocRange)/2),][order(position)]
    }
    # 只选择中心附近的 Num snp 进行分析：
    minGwasPvalueVar <- which.min(gwasEqtlInfo$pValue.gwas)
    gwasEqtlInfo <- gwasEqtlInfo[ifelse(minGwasPvalueVar>numSnp,minGwasPvalueVar-numSnp,1):ifelse( (nrow(gwasEqtlInfo)-minGwasPvalueVar)>numSnp, minGwasPvalueVar+numSnp, nrow(gwasEqtlInfo) ),]
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
  }

  if( !exists("eqtlInfo") ){
    message("Note: No eQTL associations found for gene: ", paste0(traitGenes, collapse = ", "))
    return(NULL)
  }else{
    return(list(coloc_Out_summaryAll=coloc_Out_summaryAll, gwasEqtlInfo=gwasEqtlInfo))
  }
}
