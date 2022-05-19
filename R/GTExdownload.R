
#' @title xQTLdownload_exp
#' @description
#'  download normalized gene expression at the sample level in a specified tissue.
#' @param genes Following gene types are supported:
#' \itemize{
#'   \item \strong{Gene symbol}.
#'
#'   A character string or a character vector (case ignored). like: "tp53","naDK","SDF4".
#'   \item \strong{Gencode/ensemble id} (versioned or unversioned).
#' }
#'
#' @param geneType A character string. Types of queried genes. Options: "auto" (default), "geneSymbol" and "gencodeId";
#' @param tissueSiteDetail All tissues' name can be listed with \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#'
#' @param datasetId A character string. Options: "gtex_v8" (default), "gtex_v7".
#' @param toSummarizedExperiment whether to return a data.frame or a summarizedExperiment object. Default: TRUE, return a toSummarizedExperiment object.
#' @param recordPerChunk A integer value (1-2000). number of records fetched per request (default: 150).
#' @param pathologyNotesCategories Default: pathologyNotes info is ignored.
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import stats
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom   IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @return return a SummarizedExperiment or data.frame object harbing gene expression profiles.
#' @export
#' @examples
#' \donttest{
#'   # Download gene expression with a genecode ID:
#'   expProfiles <- xQTLdownload_exp(c("ENSG00000210195.2"),
#'                                  tissueSiteDetail="Liver")
#'   # extract expression profile from SummarizedExperiment object:
#'   SummarizedExperiment::assay(expProfiles)
#'   # extract samples' detail from SummarizedExperiment object:
#'   SummarizedExperiment::colData(expProfiles)
#'
#'   # Download gene expression profiles of multiple genes:
#'   expProfiles <- xQTLdownload_exp(c("tp53","naDK","SDF4"),
#'                                  tissueSiteDetail="Artery - Coronary",
#'                                  pathologyNotesCategories=TRUE,
#'                                  toSummarizedExperiment=FALSE)
#'   expProfiles <- xQTLdownload_exp(c("tp53","naDK","SDF4"),
#'                                  tissueSiteDetail="Artery - Coronary",
#'                                  datasetId="gtex_v7")
#'
#'   # Get proteing-coding genes' expression:
#'   proT <- xQTLquery_gene (genes="protein coding")
#'   proTexp <- xQTLdownload_exp(proT$geneSymbol[1:100], tissueSiteDetail="Lung",
#'                               toSummarizedExperiment=FALSE)
#'   }
xQTLdownload_exp <- function(genes="", geneType="auto", tissueSiteDetail="Liver", datasetId="gtex_v8", toSummarizedExperiment=TRUE, recordPerChunk=150, pathologyNotesCategories=FALSE  ){
  gencodeId <- cutF <- genesUpper <- geneSymbol <- entrezGeneId <- tss <- description <- NULL
  .<-NULL
  cutNum <- recordPerChunk
  # genes=c("ENSG00000210195.2","ENSG00000078808")
  # geneType="gencodeId"
  # tissueSiteDetail="Liver"
  # datasetId="gtex_v8"
  # toSummarizedExperiment=TRUE

  # genes = miRNA$gencodeId[1000:1800]
  # geneType="gencodeId"
  # tissueSiteDetail="Liver"
  # datasetId="gtex_v8"
  # toSummarizedExperiment=TRUE
  # recordPerChunk=150

  ########## parameter check: tissueSiteDetail
  if(is.null(datasetId) ||  any(is.na(datasetId)) ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(  !(datasetId %in% c("gtex_v8", "gtex_v7"))  ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }

  ########## parameter check: tissueSiteDetail
  # import tissueSiteDetailGTEx according to datasetId
  gencodeVersion <- "v26"
  genomeBuild <- "GRCh38/hg38"
  if(datasetId == "gtex_v8"){
    # data(tissueSiteDetailGTExv8)
    # data(gencodeGeneInfoV26)
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
    # gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV26)
    gencodeVersion <- "v26"
    genomeBuild <- "GRCh38/hg38"
  }else if(datasetId == "gtex_v7"){
    # data(tissueSiteDetailGTExv7)
    # data(gencodeGeneInfoV19)
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
    # gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV19)
    gencodeVersion <- "v19"
    genomeBuild <- "GRCh37/hg19"
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be null!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }
  # convert tissueSiteDetail to tissueSiteDetailId:
  tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId


  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else if( length(genes)==1 ){
      if( genes %in% gencodeGenetype$V26 | genes %in% gencodeGenetype$V19 ){
        geneType <- "geneCategory"
      }else{
        geneType <- "geneSymbol"
      }
    }else{
      geneType <- "geneSymbol"
    }
  }


  ############ convert genes. parameter check is unnecessary for this, because xQTLquery_gene check it internally.
  message("== Validate gene :")
  geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, recordPerChunk=recordPerChunk)
  if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
    stop("gene information is null.")
  }
  message("== Done.")
  # geneInfo <-  xQTLquery_gene(c("tp53","naDK","SDF4"), "geneSymbol", "v26", "GRCh38/hg38")

  ############ get sample info:
  message("== Validate sample:")
  sampleInfo <- xQTLquery_sampleByTissue(tissueSiteDetail=tissueSiteDetail, dataType="RNASEQ", datasetId=datasetId, recordPerChunk=recordPerChunk,pathologyNotesCategories=pathologyNotesCategories )
  message("== Done.")
  if( !exists("sampleInfo") ||is.null(sampleInfo) ){
    stop("Failed to fetch sample information.")
  }

  ############### 分批下载：
  geneInfo <- cbind(geneInfo, data.table(cutF = cut(1:nrow(geneInfo),breaks=seq(0,nrow(geneInfo)+cutNum,cutNum) ) ))
  setDT(geneInfo)
  geneInfo$genesUpper <- toupper(geneInfo$genes)
  genesURL <- geneInfo[,.(genesURL=paste0(gencodeId[!is.na(gencodeId)],collapse = "%2C")),by=c("cutF")]
  if( any(unlist(lapply(genesURL$genesURL, nchar)) >3900) ){
    stop("Too many queried genes, please lower the value of \"recordPerChunk\", or reduce your input genes.")
  }
  outInfo <- data.table::data.table()
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")
  for(i in 1:nrow(genesURL)){
    # construct url:
    url1 <- paste0("https://gtexportal.org/rest/v1/expression/geneExpression?",
                   "datasetId=", datasetId,"&",
                   "gencodeId=", genesURL[i,]$genesURL, "&",
                   "tissueSiteDetailId=", tissueSiteDetailId,"&",
                   "&format=json"
    )
    url1 <- utils::URLencode(url1)
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

    tmp <- data.table::as.data.table(url1GetText2Json$geneExpression)
    if(nrow(tmp)==0){
      message("No expression profiles were found in ",tissueSiteDetail, " of thess ", length(genes), " genes!")
      message("== Done.")
      return(data.table::data.table())
    }
    tmp$genesUpper <- toupper(unlist(tmp[,geneType,with=FALSE]))
    # keep raw genes:
    tmp <- merge(geneInfo[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp, by="genesUpper",all.x=TRUE, sort = FALSE)[,-c("genesUpper")]
    outInfo <- rbind(outInfo, tmp)
    message("Downloaded part ", i, "/",nrow(genesURL),"; ", nrow(outInfo), " records.")
    rm(url1, url1GetText2Json, tmp)
  }
  # expression, fill with NA if NULL:
  tmpExp <- outInfo$data
  # tmpExp <- lapply(tmpExp, function(x){ if(is.null(x)){return( data.table::data.table( matrix(rep(NA,nrow(sampleInfo)),nrow = 1) ) )}else{return( data.table::data.table( matrix(x, nrow=1) ) )} })
  # expDT <- data.table::rbindlist(tmpExp)
  tmpExp <- lapply(tmpExp, function(x){  if(is.null(x)){ rep(NA, nrow(sampleInfo) )}else{ x } })
  expDT <- data.table::as.data.table(t(data.table::as.data.table(tmpExp)))
  data.table::setDF(expDT)
  colnames(expDT) <- sampleInfo$sampleId
  rownames(expDT) <- geneInfo$genes
  rm(tmpExp)
  geneInfo$genesUpper <- NULL
  geneInfo$cutF <- NULL
  # convert to SummarizedExperiment or not:
  if( toSummarizedExperiment ){      # construct summarizedExperiment object::
    # gene info:
    expRowRanges <- GenomicRanges::GRanges(ifelse(is.na(geneInfo$chromosome),"NA",geneInfo$chromosome),
                                           IRanges::IRanges(ifelse(is.na(geneInfo$start),-1,geneInfo$start), ifelse(is.na(geneInfo$end),-1,geneInfo$end)),
                                           strand= ifelse(is.na(geneInfo$strand),"*",geneInfo$strand),
                                           geneInfo[,.(genes, geneSymbol, gencodeId, entrezGeneId, geneType, tss, description)]
                                           # ,seqinfo = GenomeInfoDb::Seqinfo(genome= stringr::str_split(unique(includedGenes$genomeBuild)[1],stringr::fixed("/"))[[1]][2])
    )
    # sample info:
    expColData <- data.table::copy(sampleInfo)
    # meta data:
    expMetadata <- paste0("Queried ", length(genes), " genes, finally ", length(outInfo$data[[which(unlist(lapply(outInfo$data,length))>0)[1]]])," expression profiles of ",length(na.omit(outInfo$gencodeId)), " genes in ",unique(na.omit(outInfo$tissueSiteDetailId)), " were obtained. Unit: ",paste0(unique(na.omit(outInfo$unit)), collapse = ","),".")

    # outInfoSE:
    outInfoSE <- SummarizedExperiment::SummarizedExperiment(assays=list(exp=expDT),
                                                            rowRanges=expRowRanges,
                                                            colData=expColData,
                                                            metadata=expMetadata)
    message("== Done.")
    return(outInfoSE)
  }else{
    outInfoDF <- cbind(outInfo[,-c("data")], expDT)
    data.table::setDF(outInfoDF)
    message("== Done.")
    return(outInfoDF)
  }
}


#' @title xQTLdownload_eqtlSig
#' @description download significant eQTL associations of a tissue or across all tissues.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param genes A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail A character string. tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @return A data.table
#' @export
#'
#' @examples
#' \donttest{
#'  # Download eQTL info for a variant:
#'  xQTLdownload_eqtlSig("rs201327123")
#'  xQTLdownload_eqtlSig("chr1_14677_G_A_b38")
#'  xQTLdownload_eqtlSig("11_66328719_T_C_b37", datasetId="gtex_v7")
#'  xQTLdownload_eqtlSig("11_66328719_T_C_b37", datasetId="gtex_v7",
#'                        tissueSiteDetail="Skin - Sun Exposed (Lower leg)")
#'
#'  # Download eQTL association according to all tissues with genome location:
#'  varInfo <-  xQTLquery_varPos(chrom="chr1", pos=c(1102708),"gtex_v8")
#'  xQTLdownload_eqtlSig(variantName=varInfo$snpId)
#'
#'  # Download eQTL info for gene:
#'  xQTLdownload_eqtlSig(genes="ATAD3B")
#'  xQTLdownload_eqtlSig(genes=c("TP53", "SLC35E2B"), tissueSiteDetail= "Brain - Cerebellum")
#'  xQTLdownload_eqtlSig(genes="ENSG00000141510.16", datasetId="gtex_v8")
#'  xQTLdownload_eqtlSig(genes="ENSG00000141510.11", datasetId="gtex_v7",
#'                       tissueSiteDetail="Thyroid" )
#'  xQTLdownload_eqtl(genes="ENSG00000141510.11",
#'                    datasetId="gtex_v7",tissueSiteDetail="Thyroid" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  xQTLdownload_eqtlSig(variantName="rs1641513", genes="TP53", datasetId="gtex_v8")
#'  xQTLdownload_eqtlSig(variantName="rs1641513", genes="TP53", datasetId="gtex_v7")
#'  xQTLdownload_eqtlSig(variantName="chr1_1667948_A_G_b38",
#'                    genes="SLC35E2B", tissueSiteDetail="Kidney - Cortex")
#' }
xQTLdownload_eqtlSig <- function(variantName="", genes="", variantType="auto", geneType="auto", tissueSiteDetail="", datasetId="gtex_v8"){
  variantId <- snpId <- gencodeId <- geneSymbol <- pValue <- nes <- NULL
  .<-NULL
  # variantName="chr1_14677_G_A_b38"
  # variantType="variantId"
  # datasetId="gtex_v8"
  # tissueSiteDetail="Liver"
  # gene = "ENSG00000228463.9"
  # geneType="gencodeId"
  # check version:
  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if( datasetId == "gtex_v7" ){
    gencodeVersion <- "v19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }


  ##################### fetch varInfo
  if(variantName!=""){

    # auto pick variantType
    if(variantType=="auto"){
      if(stringr::str_detect(variantName, stringr::regex("^rs"))){
        variantType <- "snpId"
      }else if(stringr::str_count(variantName,"_")==4){
        variantType <- "variantId"
      }else{
        stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
      }
    }

    message("== Validate variant:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      message("The variant [",variantName, "] is not incuded in [",datasetId,"].")
      return(NULL)
    }else{
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if( all(genes !="") ){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(genes, function(g){ stringr::str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Validate gene:")
    geneInfo <- xQTLquery_gene(genes=genes, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }
  # both are null: error!
  if(variantName=="" && all(genes == "")){
    stop("Parameter \"variantName\" and \"gene\" can not be null at the same time. ")
  }

  # if( variantName is null || geneName is null; tissue 随意){
    # return all sig
  # }else if ( both is not null && tissue is not null){
    # return unsig or sig with exp
  # }
  message("== Start downloading significant eQTL associations:")
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/singleTissueEqtl?format=json",
                 "&datasetId=",datasetId,
                 ifelse(variantName=="","",paste0("&variantId=",varInfo$variantId)),
                 ifelse(all(genes==""),"",paste0("&gencodeId=",paste0(geneInfo$gencodeId, collapse = ","))),
                 ifelse(tissueSiteDetail=="","",paste0("&tissueSiteDetailId=",tissueSiteDetailId))
                 )
  url1 <- utils::URLencode(url1)
  # check network:
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  tmp <- data.table::as.data.table(url1GetText2Json$singleTissueEqtl)
  if( !exists("tmp")||nrow(tmp)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & genes!="","-",""),ifelse(genes=="","",paste0(" gene [", genes,"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail)), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  tmp <- merge(tmp, tissueSiteDetailGTEx, by = "tissueSiteDetailId")
  outInfo <- tmp[,.(variantId, snpId, gencodeId, geneSymbol, tissueSiteDetail, pValue, nes,datasetId)]
  message("=================================")
  message("Totally ", nrow(outInfo), " associatons were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & genes!=""," -",""),ifelse(all(genes==""),"",paste0(" gene: [",paste0(genes, collapse = ", "),"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail))," in ",datasetId,"."  )

  return(outInfo)
}


#' @title xQTLdownload_eqtl
#' @description download significant or unsignificant eQTL data of a tissue or across all tissues.
#'  can be quried with a gene/variant-gene pair.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned). Can not be null.
#' @param variantType A character string. "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail A character string. tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param recordPerChunk A integer value (1-500). number of records fetched per request (default: 100).
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @return A data.table
#' @export
#'
#' @examples
#' \donttest{
#'
#'  # Download eQTL info for a gene:
#'  eqtlInfo <- xQTLdownload_eqtl(gene="TP53")
#'  xQTLdownload_eqtl(gene="ATAD3B", datasetId="gtex_v7")
#'
#'  # Unversioned gencode ID in GTEx V7:
#'  eqtl_v7 <- xQTLdownload_eqtl(gene="ENSG00000141510",
#'                               datasetId="gtex_v7")
#'  # Unversioned gencode ID in GTEx V8:
#'  eqtl_v8 <- xQTLdownload_eqtl(gene="ENSG00000141510",
#'                               datasetId="gtex_v8")
#'
#'  # In a specific tissue:
#'  xQTLdownload_eqtl(gene="ENSG00000141510.11", geneType="gencodeId",
#'                       datasetId="gtex_v7", tissueSiteDetail="Thyroid" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  xQTLdownload_eqtl(variantName="rs1641513",gene="TP53", datasetId="gtex_v8")
#'  xQTLdownload_eqtl(variantName="rs11657498",gene="TP53",
#'                       datasetId="gtex_v7")
#'  xQTLdownload_eqtlSig(variantName="chr1_1667948_A_G_b38", gene="SLC35E2B",
#'                       tissueSiteDetail="Kidney - Cortex")
#'  xQTLdownload_eqtl(variantName="17_7492388_G_A_b37",gene="TP53",
#'                       tissueSiteDetail="Uterus",  datasetId="gtex_v7")
#' }
xQTLdownload_eqtl <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail="", datasetId="gtex_v8", recordPerChunk=100){
  pos <- variantId <- snpId <- gencodeId <- geneSymbol <- pValue <- nes <- se <- NULL
  .<-NULL
  querySnpId=TRUE
  # variantName="chr1_14677_G_A_b38"
  # variantType="variantId"
  # datasetId="gtex_v8"
  # tissueSiteDetail="Liver"
  # gene = "ENSG00000228463.9"
  # geneType="gencodeId"

  # variantName=""
  # variantType="variantId"
  # datasetId="gtex_v7"
  # tissueSiteDetail="Liver"
  # gene = "ATAD3B"
  # geneType="geneSymbol"

  # check version:
  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if( datasetId == "gtex_v7" ){
    gencodeVersion <- "v19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }
  # check gene:
  if( length(gene) >1 || !is.character(gene) ){
    stop("Parameter \"gene\" must be a character string. ")
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }
  # check network:
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")

  ##################### fetch geneInfo:
  if(gene !=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else if( length(gene)==1 ){
        if( gene %in% gencodeGenetype$V26 | gene %in% gencodeGenetype$V19 ){
          geneType <- "geneCategory"
        }else{
          geneType <- "geneSymbol"
        }
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Validate gene:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or leave \"gene\" as null")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"gene\" can not be null! ")
  }
  ##################### fetch varInfo
  if(variantName!=""){

    # auto pick variantType
    if(variantType=="auto"){
      if(stringr::str_detect(variantName, stringr::regex("^rs"))){
        variantType <- "snpId"
      }else if(stringr::str_count(variantName,"_")==4){
        variantType <- "variantId"
      }else{
        stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
      }
    }

    message("== Validate variant:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input, or leave \"variantName\" as null.")
    }else{
      message("== Done.")
    }
  }

  # url1 <- "https://gtexportal.org/rest/v1/association/metasoft?gencodeId=ENSG00000141510.11&datasetId=gtex_v7"
  message("== Querying significant eQTL associations from API server:")
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/metasoft?",
                 "&datasetId=",datasetId,
                 ifelse(variantName=="","",paste0("&variantId=",varInfo$variantId)),
                 "&gencodeId=",geneInfo$gencodeId
  )
  url1 <- utils::URLencode(url1)
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

  if( length(url1GetText2Json$metasoft)==0 ){
    message("No eQTL association found for gene [",gene,"]",ifelse(variantName=="","",paste0(" - variant [",variantName,"].")))
    return( data.table::data.table() )
  }else{
    message("== Done.")
  }
  tmp <- url1GetText2Json$metasoft

  tmp_info <- data.table::data.table( gencodeId=tmp$gencodeId, variantId=tmp$variantId, datasetId=tmp$datasetId, metaP=tmp$metaP )
  tmp <- tmp[,-which(names(tmp) %in% c("datasetId", "gencodeId", "metaP","variantId"))]
  outInfo <- data.table::data.table()
  for(i in 1:length(names(tmp))){
    tmp_i <- cbind( tmp_info, tmp[[i]])
    tmp_i_tissue <- tissueSiteDetailGTEx[names(tmp)[i], on="tissueSiteDetailId"]$tissueSiteDetail
    tmp_i$tissueSiteDetail <- tmp_i_tissue
    outInfo <- rbind(outInfo, tmp_i )
    outInfo <- cbind(outInfo[,-c("datasetId")], outInfo[,c("datasetId")])
  }

  if(tissueSiteDetail!=""){
    outInfo <- outInfo[tissueSiteDetail, on="tissueSiteDetail"]
  }
  # add geneSymbol:
  outInfo$geneSymbol <- geneInfo$geneSymbol
  if(querySnpId){
    # Fetch dbSNP ID:
    message("== Querying dbsnp ID from API server:")
    var_tmp <- data.table::data.table(variantId =unique(outInfo$variantId))
    var_tmp <- cbind(var_tmp, data.table::rbindlist(lapply(unique(outInfo$variantId), function(x){ splitOut = stringr::str_split(x,stringr::fixed("_"))[[1]];data.table::data.table(chrom=splitOut[1], pos=splitOut[2]) })))
    var_tmpGot <- xQTLquery_varPos(unique(var_tmp$chrom), pos=as.integer(unlist(var_tmp$pos)), datasetId = datasetId, recordPerChunk=recordPerChunk)
    outInfo <- merge(outInfo, var_tmpGot, by="variantId")
    message("== Done.")
  }else{
    outInfo$snpId <- ""
  }
  outInfo <- outInfo[,.(variantId, snpId, gencodeId, geneSymbol, tissueSiteDetail, pValue, nes, se, mValue ,datasetId)]

  message("=================================")
  message("In total of ", nrow(outInfo) ," eQTL associations were found for ", ifelse(gene=="","",paste0("Gene [", gene,"]")), ifelse(variantName=="","",paste0(" - Variant [",variantName,"]")), ifelse(tissueSiteDetail=="", paste0(" in ",length(unique(outInfo$tissueSiteDetail)),ifelse(length(unique(outInfo$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail))," in ",datasetId,"."  )
  return(outInfo)
}

#' @title This function return the gene-variant association for any given pair of gene and variant, which may or may not be significant
#' @description Only GTEx V8 are supported.
#' @param genes GENCODE ID (versioned or unversioned).
#' @param variants GTEx variant ID or dbSNP ID. Note: if variant id are provided, genome version of hg38 are required.
#' @param tissueSiteDetail  tissue site detail.
#' @param recordPerChunk A integer value (1-200). number of records fetched per request (default: 50).
#'
#' @return A data.table object
#'
#' @examples
#' \donttest{
#'    geneList = c("ENSG00000228463","ENSG00000228463.9", "ENSG00000065613.13")
#'    variants = c("rs201327123","chr1_14677_G_A_b38", "chr11_66561248_T_C_b38")
#'    xQTLdownload_eqtlPost(geneList, variants, tissueSiteDetail="Colon - Sigmoid")
#' }
# xQTLdownload_eqtlPost <- function(geneList, variantlist, tissueSiteDetail="", recordPerChunk=50){
#   # tissue convert:
#   tissueSiteDetailId <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
#   # gene split:
#   geneList <- unlist(lapply(geneList, function(x){stringr::str_split(x, fixed("."))[[1]][1]} ))
#   # construct query body:
#   queryBody <- data.frame(
#     gencodeId = geneList,
#     variantId= variantlist,
#     tissueSiteDetailId=tissueSiteDetailId )
#   queryBody$cutF <- as.character(cut(1:nrow(queryBody), breaks=seq(0, nrow(queryBody)+recordPerChunk, recordPerChunk)))
#   data.table::setDT(queryBody)
#   # query:
#   cutF <- unique(queryBody$cutF)
#   resultDT <- data.table()
#   for( i in 1:length(cutF)){
#     cutFTmp<- cutF[i]
#     message("  - Downloading eQTL asso: ", paste0(round(i/length(cutF)*100,2), "%..."))
#     queryBodyTmp <- queryBody[cutF==cutFTmp,][,c("gencodeId","variantId","tissueSiteDetailId")]
#     dataTmp <- httr::POST("https://gtexportal.org/rest/v1/association/dyneqtl", body=jsonlite::toJSON(queryBodyTmp), encode = "json")
#     dataTmp
#     resultTmp <- data.table::as.data.table( fromJSON(rawToChar(dataTmp$content))$result )
#     resultDT <- rbind(resultDT, resultTmp )
#     rm( queryBodyTmp, dataTmp,resultTmp )
#     Sys.sleep(0.1)
#   }
#
#   return(resultDT)
# }

#' @title xQTLdownload_eqtlAllAsso
#' @description download all tested variant-gene associations.
#'  source of all eQTL associations is EBI eQTL category.
#'
#' @param gene gene A gene symbol or a gencode id (versioned).
#' @param geneType geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail tissueSiteDetail A character string. tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param recordPerChunk A integer value (1-200). number of records fetched per request (default: 200).
#' @param study "gtex_v8" only GTEx v8 is supported.
#' @param withB37VariantId Whether to return the genome location(GTEx v7) of variants. Default: TRUE.
#' @import data.table
#' @import stringr
#' @return a data.table
#' @export
#'
#' @examples
#' \donttest{
#'   geneAsso <- xQTLdownload_eqtlAllAsso("ATP11B", tissueSiteDetail="Muscle - Skeletal", withB37VariantId=FALSE)
#' }
xQTLdownload_eqtlAllAsso <- function(gene="", geneType="auto", tissueSiteDetail="", recordPerChunk=250, study="gtex_v8", withB37VariantId=TRUE){
  .<-NULL
  variantId <- variant <- b37VariantId <- snpId <- NULL
  # gene="CYP2W1"
  # geneType="geneSymbol"
  # tissueSiteDetail="Lung"

  if(toupper(study)=="GTEX_V8"){
    study="GTEx_V8"
  }

  ########## check genes
  if( is.null(gene) ||  any(is.na(gene)) || any(gene=="") ||length(gene)==0 ){
    stop("Parameter \"gene\" can not be NULL or NA!")
  }else if( length(gene)!=1 ){
    stop("Gene list is not supported, please input a gene symbol or gencode ID.")
  }

  ########## geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("auto","geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"auto\", \"geneSymbol\", and \"gencodeId\".")
  }

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  ########## parameter check: tissueSiteDetail
  # import tissueSiteDetailGTEx according to datasetId
  tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  gencodeVersion <- "v26"
  qtl_groups <- EBIquery_allTerm("qtl_groups")
  qtl_tissue <- merge( tissueSiteDetailGTEx,qtl_groups, by.x="tissueSiteDetailId", by.y="qtl_group")
  setDT(qtl_tissue)
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) || tissueSiteDetail==""   ){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if( !(tissueSiteDetail %in% c("All", tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c(paste0(1:nrow(qtl_tissue),". ",qtl_tissue$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }
  # convert tissueSiteDetail to tissueSiteDetailId:
  tissueSiteDetailId <- qtl_tissue[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId

  ##################### fetch geneInfo:
  if(gene !=""){
    message("== Validate gene:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion)
    geneInfo <- na.omit(geneInfo)
    # geneInfoV19 <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v19")
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or set gene with gencodeId.")
    }else{
      message("== Done.")
    }

    # else if( nrow(geneInfo)>1 || nrow(geneInfoV19)>1 ){
    #   stop("Totally, ",nrow(geneInfo), " gencode ID of queried gene [", gene,"] were detected, please enter the gencode ID (versioned or unversioned) for querying!")
    # }

  }else{
    stop("Parameter \"gene\" can not be null! ")
  }
  geneInfo$gencodeIdUnv <-stringr::str_split(geneInfo$gencodeId, stringr::fixed("."))[[1]][1]

  # construct url:
  url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/studies/",study,
                 "/associations?links=False&gene_id=", geneInfo$gencodeIdUnv,
                 "&qtl_group=",tissueSiteDetailId)
  # check network:
  bestFetchMethod <- apiEbi_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    message("Note: EBI API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  gtexAsoo <- fetchContentEbi(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2],  termSize=1000)
  gtexAsooList <- do.call(c, gtexAsoo)
  if(length(gtexAsooList)==0){
    message("No association found!")
    return(NULL)
  }
  gtexAsooList <- lapply(gtexAsooList, function(x){ data.table(variantId= x$variant, snpId=x$rsid,type=x$type,maf=x$maf,beta=x$beta,
                                                               chrom=x$chromosome, pos=x$position,
                                                               # ref=x$ref, alt=x$alt, gene_id= x$gene_id, study_id=x$study_id
                                                               se=x$se, median_tpm=x$median_tpm, pValue=x$pvalue, totalAlleles=x$an, allelCounts=x$ac, imputationR2=x$r2,
                                                               tissue_label=x$tissue_label,tissue=x$tissue,qtl_group=x$qtl_group, condition=x$condition,
                                                               molecular_trait_id = x$molecular_trait_id
                                                               ) })
  gtexAsooDT <- data.table::rbindlist(gtexAsooList, fill=TRUE)
  gtexAsooDT$geneSymbol <- geneInfo$geneSymbol
  gtexAsooDT$gencodeId_GTEX_v8 <- geneInfo$gencodeId
  gtexAsooDT$gencodeId_GTEX_v7 <- ifelse(nrow(geneInfoV19)>0, geneInfoV19$gencodeId, "")
  if(withB37VariantId){
    # add dbSNP id and  hg19 cordinate:
    gtexAsooDTb37 <- xQTLquery_varPos(chrom = paste0("chr",unique(gtexAsooDT$chrom)), pos = gtexAsooDT$pos, datasetId = "gtex_v8", recordPerChunk = recordPerChunk)
    gtexAsooDTb37$variantId <- unlist(lapply(gtexAsooDTb37$variantId, function(x){ splitInfo=stringr::str_split(x, stringr::fixed("_"))[[1]]; paste0(splitInfo[-5], collapse="_") }))
    gtexAsooDT <- merge(gtexAsooDT, gtexAsooDTb37[,.(variantId, b37VariantId)], by=c("variantId"), all.x=TRUE )
    gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
    gtexAsooDT <- cbind(gtexAsooDT[,.(snpId, variantId, b37VariantId)], gtexAsooDT[,-c("snpId", "variantId", "b37VariantId", "chrom", "pos")])
    return(gtexAsooDT)
  }else{
    return(gtexAsooDT)
  }
}

#' @title Fetch all eQTL associations from EBI eQTL category.
#'
#' @param gene gene symbol or gencode ID.
#' @param geneType A character string. "geneSymbol"(default), "gencodeId" or "geneCategory".
#' @param refSeq A character string.
#' @param rsid dbSNP ID.
#' @param variantId A character string. like "chr17_5294580_CT_C_b38".
#' @param pThreshold Float. set a threshold for p-value to filter the result. Default: 1, aQTL associations that have a p-value<1 will be returned.
#' @param tissueSiteDetail tissue.
#'
#' @examples
#' \donttest{
#'   xQTLdownload_aQTLAllAsso( gene="TP53", tissueSiteDetail = "Adipose - Visceral (Omentum)" )
#'   xQTLdownload_aQTLAllAsso( refSeq = "NM_001291581.2", tissueSiteDetail = "Adipose - Visceral (Omentum)" )
#'   xQTLdownload_aQTLAllAsso( rsid = "rs3026133", tissueSiteDetail = "Adipose - Visceral (Omentum)", pThreshold=1e-1 )
#'   xQTLdownload_aQTLAllAsso( variantId = "chr17_5289661_A_G_b38", tissueSiteDetail = "Adipose - Visceral (Omentum)" )
#' }
# xQTLdownload_aQTLAllAsso<- function(gene="", geneType="geneSymbol", refSeq="", rsid="",  variantId="", pThreshold=1, tissueSiteDetail=""){
  # gene="RABEP1"
  # geneType="geneSymbol"
  # refSeq = "NM_001291581.2"
  # tissueSiteDetail <- "Adipose - Visceral (Omentum)"


  # if(tissueSiteDetail == ""){
  #   stop( "Tissue ",tissueSiteDetail, " not found!"  )
  # }
  # tissueSiteDetailId  <- tissueSiteDetailGTExv8[tissueSiteDetail, on="tissueSiteDetail"]$tissueSiteDetailId
  # if(gene !=""){
  #   geneInfo <- xQTLquery_gene(gene, geneType = geneType)
  #   # create sql:
  #   sqlForQuery <- paste0("select * from ", paste0("aQTL_",tissueSiteDetailId), " where geneSymbol = \'",geneInfo$geneSymbol,"\' AND `p-value`<",pThreshold )
  #   con <- DBI::dbConnect(RMySQL::MySQL(), dbname = 'GTExDB', host = "172.18.200.246", port = 3306, user = "GTExAnonymous", password = "GTEx_Anonymous123")
  #   # # extract data:
  #   message("== mysql connected, downloading... ")
  #   aqltInfo <- dbGetQuery(con, sqlForQuery)
  #   dbDisconnect(con)
  #   message("== done.")
  #   # message("== mysql closed || ", date(), " || [", nrow(snpLD), "] LD infor obtained!")
  #   return(aqltInfo)
  # }else if(refSeq!=""){
  #   sqlForQuery <- paste0("select * from ", paste0("aQTL_",tissueSiteDetailId), " where refSeq = \'",refSeq,"\' AND `p-value`<",pThreshold )
  #   con <- DBI::dbConnect(RMySQL::MySQL(), dbname = 'GTExDB', host = "172.18.200.246", port = 3306, user = "GTExAnonymous", password = "GTEx_Anonymous123")
  #   message("== mysql connected, downloading... ")
  #   aqltInfo <- dbGetQuery(con, sqlForQuery)
  #   dbDisconnect(con)
  #   message("== done.")
  #   return(aqltInfo)
  # }else if(rsid!=""){
  #   sqlForQuery <- paste0("select * from ", paste0("aQTL_",tissueSiteDetailId), " where rsid = \'",rsid,"\' AND `p-value`<",pThreshold )
  #   con <- DBI::dbConnect(RMySQL::MySQL(), dbname = 'GTExDB', host = "172.18.200.246", port = 3306, user = "GTExAnonymous", password = "GTEx_Anonymous123")
  #   message("== mysql connected, downloading... ")
  #   aqltInfo <- dbGetQuery(con, sqlForQuery)
  #   dbDisconnect(con)
  #   message("== done.")
  #   return(aqltInfo)
  # }else if(variantId!=""){
  #   sqlForQuery <- paste0("select * from ", paste0("aQTL_",tissueSiteDetailId), " where SNP = \'",variantId,"\' AND `p-value`<",pThreshold )
  #   con <- DBI::dbConnect(RMySQL::MySQL(), dbname = 'GTExDB', host = "172.18.200.246", port = 3306, user = "GTExAnonymous", password = "GTEx_Anonymous123")
  #   message("== mysql connected, downloading... ")
  #   aqltInfo <- dbGetQuery(con, sqlForQuery)
  #   dbDisconnect(con)
  #   message("== done.")
  #   return(aqltInfo)
  # }else{
  #   stop("Please check your input!")
  # }

# }


#' @title xQTLdownload_sqtlSig
#' @description download significant sQTL associations of a tissue or across all tissues.
#'  Only GTEx v8 is supported.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail A character string. tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @return A data.table
#' @export
#'
#' @examples
#' \donttest{
#'  # Download sQTL info for a variant:
#'  xQTLdownload_sqtlSig(variantName="rs201327123")
#'  xQTLdownload_sqtlSig(variantName="chr1_14677_G_A_b38")
#'  xQTLdownload_sqtlSig(variantName="chr1_14677_G_A_b38",
#'                       tissueSiteDetail="Whole Blood")
#'
#'  # Download sQTL association according to all tissues with genome location:
#'  xQTLquery_varPos(chrom="chr1", pos=c(1102708),"gtex_v8")
#'  xQTLdownload_sqtlSig(variantName=varInfo$snpId, variantType="snpId")
#'
#'  # Download sQTL info for a gene:
#'  xQTLdownload_sqtlSig(gene="ATAD3B")
#'  xQTLdownload_sqtlSig(gene="ENSG00000141510.16")
#'  xQTLdownload_sqtlSig(gene="ENSG00000141510.16",
#'                       tissueSiteDetail="Lung" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  xQTLdownload_sqtlSig(variantName="rs546057177", gene="TP53")
#'  xQTLdownload_sqtlSig(variantName="chr17_7465085_A_G_b38",
#'                       gene="TP53", tissueSiteDetail="Lung")
#' }
xQTLdownload_sqtlSig <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail="" ){
  .<-NULL
  variantId <- snpId <- gencodeId <- geneSymbol <- phenotypeId <- pValue <- nes <- NULL
  # variantName="chr1_739465_TTTTG_T_b38"
  # variantType="variantId"
  # variantName = "rs1450891501"
  # variantType="snpId"
  # gene="TP53"
  # geneType="geneSymbol"
  # gene="ENSG00000141510.16"
  # geneType = "gencodeId"
  # tissueSiteDetail=""
  # datasetId= "gtex_v8"

  datasetId="gtex_v8"
  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else{
    stop("sQTL only support \"gtex_v8\" of datasetId.")
  }
  # check gene:
  if( length(gene) >1 || !is.character(gene) ){
    stop("Parameter \"gene\" must be a character string. ")
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }

  ##################### fetch varInfo
  if( variantName!=""){

    # auto pick variantType
    if(variantType=="auto"){
      if(stringr::str_detect(variantName, stringr::regex("^rs"))){
        variantType <- "snpId"
      }else if(stringr::str_count(variantName,"_")==4){
        variantType <- "variantId"
      }else{
        stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
      }
    }

    message("== Validate variant:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if(gene !=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Validate gene:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or leave \"gene\" as null")
    }else{
      message("== Done.")
    }
  }
  # both are null: error!
  if(variantName=="" && gene == ""){
    stop("Parameter \"variantName\" and \"gene\" can not be null at the same time. ")
  }

  # if( variantName is null || geneName is null; tissue 随意){
  # return all sig
  # }else if ( both is not null && tissue is not null){
  # return unsig or sig with exp
  # }
  message("== Querying significant sQTL associations from API server:")
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/singleTissueSqtl?format=json",
                 "&datasetId=",datasetId,
                 ifelse(variantName=="","",paste0("&variantId=",varInfo$variantId)),
                 ifelse(gene=="","",paste0("&gencodeId=",geneInfo$gencodeId)),
                 ifelse(tissueSiteDetail=="","",paste0("&tissueSiteDetailId=",tissueSiteDetailId))
  )
  # check network:
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  url1 <- utils::URLencode(url1)
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  tmp <- data.table::as.data.table(url1GetText2Json$singleTissueSqtl)
  if(nrow(tmp)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & gene!="","-",""),ifelse(gene=="","",paste0(" gene [", gene,"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail)), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  tmp <- merge(tmp, tissueSiteDetailGTEx, by = "tissueSiteDetailId")
  outInfo <- tmp[,.(variantId, snpId, gencodeId, geneSymbol, phenotypeId, tissueSiteDetail, pValue, nes, datasetId)]
  message("=================================")
  message("Totally ", nrow(outInfo), " associatons were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & gene!=""," -",""),ifelse(gene=="","",paste0(" gene [", gene,"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail))," in ",datasetId,"."  )

  return(outInfo)
}

#' @title xQTLdownload_eqtlExp
#' @description download normalized expression of gene for a eQTL pair.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail A character string. tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import stats
#' @import tidyr
#' @return A data.table
#' @export
#' @examples
#' \donttest{
#'  # Download exp in different tissues:
#'  xQTLdownload_eqtlExp(variantName="rs1641513",gene="TP53", tissueSiteDetail="Liver")
#'  xQTLdownload_eqtlExp(variantName="rs1641513",gene="ATAD3B",
#'                       tissueSiteDetail="Lung", datasetId="gtex_v8")
#'
#'  # Download exp in gtex v7:
#'  xQTLdownload_eqtlExp(variantName="rs140894808",gene="ATAD3B",
#'                       tissueSiteDetail="Adipose - Visceral (Omentum)",
#'                       datasetId="gtex_v7")
#'
#'  # Dowload exp using variant ID and gencode ID.
#'  xQTLdownload_eqtlExp(variantName="chr1_14677_G_A_b38",gene="ENSG00000228463.9",
#'                      tissueSiteDetail="Stomach")
#' }
xQTLdownload_eqtlExp <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail="", datasetId="gtex_v8"){
  # variantName="chr1_14677_G_A_b38"
  # variantType="variantId"
  # datasetId="gtex_v8"
  # tissueSiteDetail="Liver"
  # gene = "ENSG00000228463.9"
  # geneType="gencodeId"
  # check version:
  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if( datasetId == "gtex_v7" ){
    gencodeVersion <- "v19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }
  # check gene:
  if( length(gene) >1 || !is.character(gene) ){
    stop("Parameter \"gene\" must be a character string. ")
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if(  !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }else{
    message("Invalid parameter \"tissueSiteDetail\".")
    return(data.table::data.table())
  }

  ##################### fetch varInfo
  if(variantName!=""){

    # auto pick variantType
    if(variantType=="auto"){
      if(stringr::str_detect(variantName, stringr::regex("^rs"))){
        variantType <- "snpId"
      }else if(stringr::str_count(variantName,"_")==4){
        variantType <- "variantId"
      }else{
        stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
      }
    }

    message("== Validate variant:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo) || !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"variantName\" can not be null!")
  }

  ##################### fetch geneInfo:
  if(gene !=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Validate gene:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    if(nrow(geneInfo)==0|| is.null(geneInfo) || !exists("geneInfo")){
      stop("Invalid gene name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"gene\" can not be null!")
  }

  message("== Downloading expression...")
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # construct url
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/dyneqtl?",
                 "gencodeId=",geneInfo$gencodeId,"&",
                 "variantId=",varInfo$variantId,"&",
                 "tissueSiteDetailId=",tissueSiteDetailId,"&",
                 "datasetId=",datasetId
  )
  url1 <- utils::URLencode(url1)
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  expData <- data.table::data.table(normExp=url1GetText2Json$data,genotypes=url1GetText2Json$genotypes)

  if(nrow(expData)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & gene!="","-",""),ifelse(gene=="","",paste0(" gene [", gene,"]")), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  message("=================================")
  message("== Summary: ")
  message("[ pValue ]:          ",url1GetText2Json$pValue)
  message("[ pValueThreshold ]: ",url1GetText2Json$pValueThreshold)
  message("[ nes ]:             ",url1GetText2Json$nes)
  message("[ maf ]:             ",url1GetText2Json$maf)
  message("[ error ]:           ",url1GetText2Json$error)
  message("== Normalized expression and genotypes of [", nrow(expData), "] samples were found for", paste0(" variant [", variantName,"]"), " -",paste0(" gene [", gene,"]")," pair in tissue [", tissueSiteDetail,"] in [",datasetId,"]."  )
  message("== Genotype: ",url1GetText2Json$hetCount," het; ",url1GetText2Json$homoAltCount," hom; " ,url1GetText2Json$homoRefCount," ref.")

  message("== For more normalization method details, you can visit: https://www.gtexportal.org/home/faq#normalization. ")
  return(expData)
}


#' @title xQTLdownload_sqtlExp
#' @description
#'  download normalized expression of intron for a sQTL pair.
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param phenotypeId A character string. Format like: "chr1:497299:498399:clu_54863:ENSG00000239906.1"
#' @param variantType A character string. "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueSiteDetail A character string. Tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import stats
#' @import tidyr
#' @return A data.table
#' @export
#'
#' @examples
#' \donttest{
#'  # Download exp in different tissues:
#'  xQTLdownload_sqtlExp(variantName="rs1450891501",
#'                       phenotypeId="chr1:497299:498399:clu_54863:ENSG00000239906.1",
#'                       tissueSiteDetail="Lung")
#'
#'  # Dowload exp using variant ID.
#'  xQTLdownload_sqtlExp(variantName="chr1_1259424_T_C_b38",
#'                       phenotypeId=" chr1:1487914:1489204:clu_52051:ENSG00000160072.19",
#'                       tissueSiteDetail="Adipose - Subcutaneous")
#' }
xQTLdownload_sqtlExp <- function(variantName="", phenotypeId="", variantType="auto", tissueSiteDetail="", datasetId="gtex_v8"){
  # variantName="chr1_739465_TTTTG_T_b38"
  # phenotypeId="chr1:497299:498399:clu_54863:ENSG00000239906.1"
  # variantType="variantId"
  # tissueSiteDetail="Lung"
  # datasetId="gtex_v8"
  # check version:
  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if( datasetId == "gtex_v7" ){
    gencodeVersion <- "v19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }
  # check gene:
  if( length(phenotypeId) >1 || !is.character(phenotypeId) ){
    stop("Parameter \"phenotypeId\" must be a character string. ")
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if(  !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }else{
    message("Invalid parameter \"tissueSiteDetail\".")
    return(data.table::data.table())
  }

  ##################### fetch varInfo
  if(variantName!=""){

    # auto pick variantType
    if(variantType=="auto"){
      if(stringr::str_detect(variantName, stringr::regex("^rs"))){
        variantType <- "snpId"
      }else if(stringr::str_count(variantName,"_")==4){
        variantType <- "variantId"
      }else{
        stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
      }
    }

    message("== Validate variant:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo) || !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"variantName\" can not be null!")
  }

  ##################### fetch geneInfo:
  # if(gene !=""){
  #   message("== Querying gene info from API server:")
  #   geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
  #   if(nrow(geneInfo)==0|| is.null(geneInfo) || !exists("geneInfo")){
  #     stop("Invalid gene name or type, please correct your input.")
  #   }else{
  #     message("== Done.")
  #   }
  # }else{
  #   stop("Parameter \"gene\" can not be null!")
  # }

  message("== Downloading expression...")
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # construct url
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/dynsqtl?",
                 "phenotypeId=",phenotypeId,"&",
                 "variantId=",varInfo$variantId,"&",
                 "tissueSiteDetailId=",tissueSiteDetailId,"&",
                 "datasetId=",datasetId
  )
  url1 <- utils::URLencode(url1)
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  expData <- data.table::data.table(normExp=url1GetText2Json$data,genotypes=url1GetText2Json$genotypes)

  if(nrow(expData)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & phenotypeId!="","-",""),ifelse(phenotypeId=="","",paste0(" phenotypeId [", phenotypeId,"]")), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  message("=================================")
  message("== Summary: ")
  message("[ pValue ]:          ",url1GetText2Json$pValue)
  message("[ pValueThreshold ]: ",url1GetText2Json$pValueThreshold)
  message("[ nes ]:             ",url1GetText2Json$nes)
  message("[ maf ]:             ",url1GetText2Json$maf)
  message("[ error ]:           ",url1GetText2Json$error)
  message("== Normalized expression and genotypes of [", nrow(expData), "] samples were found for", paste0(" variant [", variantName,"]"), " -",paste0(" phenotypeId [", phenotypeId,"]")," pair in tissue [", tissueSiteDetail,"] in [",datasetId,"]."  )
  message("== Genotype: ",url1GetText2Json$hetCount," het; ",url1GetText2Json$homoAltCount," hom; " ,url1GetText2Json$homoRefCount," ref.")

  message("== For more normalization method details, you can visit: https://www.gtexportal.org/home/faq#normalization. ")
  return(expData)
}

#' @title xQTLdownload_ld
#' @description download linkage disequilibrium data of the variants associated with this gene.
#'
#' @param gene A gene symbol, gencode id (versioned), or a charater string of gene type.
#' @param geneType  A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId A character string. Options: "gtex_v8" (default), "gtex_v7".
#' @param recordPerChunk A integer value (1-500). number of records fetched per request (default: 100).
#' @import data.table
#' @import stringr
#' @import utils
#' @import curl
#' @import jsonlite
#' @import tidyr
#' @return A data.frame
#' @export
#'
#' @examples
#' \donttest{
#'  xQTLdownload_ld("TP53" )
#'  xQTLdownload_ld("TP53", datasetId="gtex_v7")
#'  xQTLdownload_ld(gene="ENSG00000008128.22")
#' }
xQTLdownload_ld <- function(gene = "", geneType="geneSymbol", datasetId = "gtex_v8", recordPerChunk=300){
  .<-NULL
  variantId<-snpId <- snpId_1 <- variantId_1 <- snpId_2 <-variantId_2<- ldScore <- NULL
  # check genes
  if( is.null(gene) ||  any(is.na(gene)) || any(gene=="") ||length(gene)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }else if(length(gene)!=1){
    stop("Parameter \"genes\" should be a character string!")
  }

  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  # parameter check: datasetId
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }

  # Fetch gene info:
  message("== Validate gene:")
  geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType, gencodeVersion=gencodeVersion)
  if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
    stop("The gene [",gene,"] you entered could not be found!")
  }
  message("== Done.")

  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  url1 <- paste0("https://gtexportal.org/rest/v1/dataset/ld?format=json&",
                 "gencodeId=",geneInfo$gencodeId,
                 "&datasetId=", datasetId)
  url1 <- utils::URLencode(url1)
  message("== Downloading ld...")
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  ldInfo <- data.table::data.table(url1GetText2Json$ld)
  if(nrow(ldInfo)==0){
    message("No LD information were found for", ifelse(gene=="","",paste0(" gene [", gene,"]")), " in ",datasetId)
    return(NULL)
  }
  colnames(ldInfo) <- c("variant1_variant2", "ldScore")
  # ldInfo <- cbind(ldInfo, rbindlist(lapply(ldInfo$variant1_variant2, function(x){ splitOut = stringr::str_split_fixed(x, stringr::fixed(","),2); data.table( variantId_1=splitOut[1], variantId_2=splitOut[2] ) })))
  ldInfo <- tidyr::separate(ldInfo, "variant1_variant2", into=c("variantId_1","variantId_2"),sep=",",)

  # Fetch SNP ID:
  variantInfo <- data.table::data.table(variantId = union(ldInfo$variantId_1, ldInfo$variantId_2))
  variantInfo <- cbind(variantInfo, tidyr::separate(variantInfo, "variantId", into=c("chrom", "pos", "ref", "alt", "genome"), sep="_"))
  variantInfo <- xQTLquery_varPos(chrom = unique(variantInfo$chrom), pos = as.integer(variantInfo$pos), datasetId = datasetId, recordPerChunk=recordPerChunk)

  # merge:
  ldInfo <- merge(ldInfo, variantInfo[,.(variantId, snpId)], by.x="variantId_1", by.y = "variantId", all.x = TRUE)
  ldInfo <- merge(ldInfo, variantInfo[,.(variantId, snpId)], by.x="variantId_2", by.y = "variantId", all.x = TRUE, suffixes = c("_1","_2"))

  ldInfo$ldScore <- as.numeric(ldInfo$ldScore)
  ldInfo <- ldInfo[,.(snpId_1, variantId_1, snpId_2, variantId_2, ldScore)][order(-ldScore,  snpId_1, snpId_2),]
  message("== Totally, ",nrow(ldInfo)," LD relations for ",nrow(variantInfo) ," variants of gene [",gene,"] were found in ", datasetId)
  return(ldInfo)
}


#' @title xQTLdownload_egene
#' @description download eGenes (eQTL Genes).
#'  eGenes are genes that have at least one significant cis-eQTL acting upon them. Results may be filtered by tissue.
#' @param gene  A charater string of gene symbol, gencode id (versioned).
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @param tissueSiteDetail A character string. Tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param recordPerChunk A integer value (1-200). number of records fetched per request (default: 200).
#' @import data.table
#' @import stringr
#' @import utils
#' @return a data.table
#' @export
#'
#' @examples
#' \donttest{
#'  eGeneInfoAlltissue <- xQTLdownload_egene()
#'  eGeneInfo <- xQTLdownload_egene("TP53")
#'  eGeneInfo <- xQTLdownload_egene(tissueSiteDetail="Lung", recordPerChunk=2000)
#'  eGeneInfo <- xQTLdownload_egene("ENSG00000141510.16")
#'  eGeneInfo <- xQTLdownload_egene("DDX11", datasetId="gtex_v7", tissueSiteDetail="Artery - Tibial" )
#' }
xQTLdownload_egene <- function(gene = "", geneType="auto", datasetId = "gtex_v8", tissueSiteDetail="", recordPerChunk=200){
  .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss <- log2AllelicFoldChange <- empiricalPValue <- pValue <- pValueThreshold <- qValue <-NULL
  # gene="DDX11"
  # geneType="geneSymbol"
  # datasetId = "gtex_v8"
  # tissueSiteDetail="Lung"
  # recordPerChunk=100

  page_tmp <- 0
  pageSize_tmp <- recordPerChunk

  # parameter check: datasetId
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail))){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }

  # Fetch gene info:
  if(gene!=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Validate gene:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType, gencodeVersion=gencodeVersion)
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",gene,"] you entered could not be found!")
    }
    message("== Done.")
  }

  message("Downloading eGenes..")
  # url1 <- "https://gtexportal.org/rest/v1/association/egene?pageSize=250&searchTerm=ENSG00000013573.16&sortBy=log2AllelicFoldChange&sortDirection=asc&tissueSiteDetailId=Thyroid&datasetId=gtex_v8"
  outInfo <- data.table::data.table()
  url1 <- paste0("https://gtexportal.org/rest/v1/association/egene?",
                 "page=",page_tmp,"&",
                 "pageSize=",pageSize_tmp,
                 ifelse(gene=="","",paste0("&searchTerm=",geneInfo$gencodeId)),
                 "&sortBy=log2AllelicFoldChange&sortDirection=asc",
                 ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                 "&datasetId=",datasetId)
  url1 <- utils::URLencode(url1)
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")
  # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  tmp <- data.table::as.data.table(url1GetText2Json$egene)
  outInfo <- rbind(outInfo, tmp)
  message("Records: ",nrow(outInfo),"/",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
  page_tmp<-page_tmp+1
  while( page_tmp <= (url1GetText2Json$numPages-1) ){
    url1 <- paste0("https://gtexportal.org/rest/v1/association/egene?",
                   "page=",page_tmp,"&",
                   "pageSize=",pageSize_tmp,
                   ifelse(gene=="","&",paste0("&searchTerm=",geneInfo$gencodeId)),
                   "sortBy=log2AllelicFoldChange&sortDirection=asc",
                   ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                   "&datasetId=",datasetId)
    url1 <- utils::URLencode(url1)
    # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    tmp <- data.table::as.data.table(url1GetText2Json$egene)
    outInfo <- rbind(outInfo, tmp)
    message("Records: ",nrow(outInfo),"/",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp <- page_tmp+1
  }
  outInfo <- merge(outInfo, tissueSiteDetailGTEx, by="tissueSiteDetailId")

  if(gene!=""){
    outInfo <- merge( geneInfo[,.(gencodeId, geneSymbol, entrezGeneId, geneType, chromosome, start, end, tss)],
                      outInfo[,.(gencodeId, log2AllelicFoldChange, empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)],
                      by="gencodeId")
  }else{
    outInfo <- outInfo[,.(gencodeId, log2AllelicFoldChange, empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)]
  }

  return(outInfo)

}


#' @title xQTLdownload_sgene
#' @description download sGenes (sQTL Genes).
#'  sGenes are genes that have at least one significant sQTL acting upon them. Results may be filtered by tissue.
#' @param gene  A charater string of gene symbol, gencode id (versioned). Can be null.
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId A character string. only support "gtex_v8". Default: "gtex_v8".
#' @param tissueSiteDetail A character string. Tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param recordPerChunk A integer value (1-2000). number of records fetched per request (default: 2000).
#' @import data.table
#' @import stringr
#' @import utils
#' @return a data.table
#' @export
#'
#' @examples
#' \donttest{
#'  # don't run:
#'  #sGeneInfoAlltissue <- xQTLdownload_sgene()
#'
#'  sGeneInfo <- xQTLdownload_sgene(tissueSiteDetail="Lung", recordPerChunk=2000)
#'  sGeneInfo <- xQTLdownload_sgene("ENSG00000141510.16",  tissueSiteDetail="Lung")
#'  sGeneInfo <- xQTLdownload_sgene("DDX11", tissueSiteDetail="Artery - Tibial" )
#' }
xQTLdownload_sgene <- function(gene = "", geneType="auto", datasetId = "gtex_v8", tissueSiteDetail="", recordPerChunk=2000){
  .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss <- log2AllelicFoldChange <- empiricalPValue <- pValue <- pValueThreshold <- qValue <-NULL
  # gene="ENSG00000013573.16"
  # geneType="gencodeId"
  # datasetId = "gtex_v8"
  # tissueSiteDetail="Thyroid"
  # recordPerChunk=2000

  page_tmp <- 0
  pageSize_tmp <- recordPerChunk


  # parameter check: datasetId
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail))){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }

  # Fetch gene info:
  if(gene!=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Validate gene:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType, gencodeVersion=gencodeVersion)
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",gene,"] you entered could not be found!")
    }
    message("== Done.")
  }

  message("Downloading sgenes...")

  # url1 <- "https://gtexportal.org/rest/v1/association/egene?pageSize=250&searchTerm=ENSG00000013573.16&sortBy=log2AllelicFoldChange&sortDirection=asc&tissueSiteDetailId=Thyroid&datasetId=gtex_v8"
  outInfo <- data.table::data.table()
  url1 <- paste0("https://gtexportal.org/rest/v1/association/sgene?",
                 "page=",page_tmp,"&",
                 "pageSize=",pageSize_tmp,
                 "&sortBy=empiricalPValue&sortDirection=asc&",
                 "datasetId=",datasetId,
                 ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                 ifelse(gene=="","",paste0("&gencodeId=",geneInfo$gencodeId))
                 )
  # message(url1)
  url1 <- utils::URLencode(url1)
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")
  # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  tmp <- data.table::as.data.table(url1GetText2Json$sgene)
  outInfo <- rbind(outInfo, tmp)
  message("Records: ",nrow(outInfo),"/",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
  page_tmp<-page_tmp+1
  while( page_tmp <= (url1GetText2Json$numPages-1) ){
    url1 <- paste0("https://gtexportal.org/rest/v1/association/sgene?",
                   "page=",page_tmp,"&",
                   "pageSize=",pageSize_tmp,
                   "&sortBy=empiricalPValue&sortDirection=asc&",
                   "datasetId=",datasetId,
                   ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                   ifelse(gene=="","",paste0("&gencodeId=",geneInfo$gencodeId))
    )
    url1 <- utils::URLencode(url1)
    # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    tmp <- data.table::as.data.table(url1GetText2Json$sgene)
    outInfo <- rbind(outInfo, tmp)
    message("Records: ",nrow(outInfo),"/",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp <- page_tmp+1
  }
  outInfo <- merge(outInfo, tissueSiteDetailGTEx, by="tissueSiteDetailId")

  if(gene!=""){
    outInfo <- merge( geneInfo[,.(gencodeId, geneSymbol, entrezGeneId, geneType, chromosome, start, end, tss)],
                      outInfo[,.(gencodeId, phenotypeId, nPhenotypes, empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)],
                      by="gencodeId")
  }else{
    outInfo <- outInfo[,.(gencodeId, phenotypeId, nPhenotypes,empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)]
  }

  return(outInfo)

}


#' @title xQTLdownload_geneMedExp
#' @description download genes' median expression in a tissue or across all tissues.
#'
#' @param genes A character vector of gene symbol/gencode ID.
#' @param geneType A character string. "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId "gtex_v7" or "gtex_v8"(default)
#' @param tissueSiteDetail A character string. Tissue detail can be listed using \"tissueSiteDetailGTExv8\" or \"tissueSiteDetailGTExv7\"
#' @param recordPerChunk 1-2000. Default: 150.
#' @import data.table
#' @import utils
#' @import stringr
#' @return A data.frame
#' @export
#'
#' @examples
#' \donttest{
#'  geneMedExp <- xQTLdownload_geneMedExp(genes="TP53")
#'  geneMedExp <- xQTLdownload_geneMedExp(genes=c("TP53", "IRF5"))
#' }
xQTLdownload_geneMedExp <- function(genes="", geneType="auto", datasetId="gtex_v8", tissueSiteDetail="", recordPerChunk=150 ){
  .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss<-strand <- NULL
  # check genes
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }

  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("auto","geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"auto\",\"geneSymbol\", \"gencodeId\".")
  }


  # parameter check: datasetId
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail))){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }

  # Fetch gene info:
  if(genes!="" || length(genes)>0){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Validate gene:")
    suppressMessages( geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, recordPerChunk = recordPerChunk))
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",stringr::str_sub(paste(genes, collapse = ","),1,20),"....] you entered could not be found!")
    }
    message("== Done.")
  }

  message("Downloading gene median expression...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))


  #
  outInfo <- data.table::data.table()
  genesCut <- data.table::data.table(gencodeId=geneInfo$gencodeId, ID=1:nrow(geneInfo), cutF = as.character(cut(1:nrow(geneInfo),breaks=seq(0,nrow(geneInfo)+recordPerChunk,recordPerChunk) )) )
  genesURL <- genesCut[,.(genesURL=paste0(gencodeId,collapse = "%2C")),by=c("cutF")]
  if( any(unlist(lapply(genesURL$genesURL, nchar)) >3900) ){
    stop("Too many queried genes, please lower the value of \"recordPerChunk\", or reduce your input genes.")
  }
  outInfo <- data.table::data.table()
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")
  for(i in 1:nrow(genesURL)){
    # construct url:
    url1 <- paste0("https://gtexportal.org/rest/v1/expression/medianGeneExpression?",
                   "datasetId=",datasetId,"&",
                   "gencodeId=", genesURL[i,]$genesURL,
                   ifelse(tissueSiteDetail=="","&", paste0("&tissueSiteDetailId=", tissueSiteDetailId)),
                   "format=json"
    )
    url1 <- utils::URLencode(url1)
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$medianGeneExpression)
    if( nrow(url1GetText2Json2DT)==0 ){
      message( "0 record fatched!" )
      return(NULL)
    }
    tmp <- url1GetText2Json2DT[,.(gencodeId, geneSymbol, median, unit,tissueSiteDetailId)]

    outInfo <- rbind(outInfo, tmp)
    message("Downloaded  ", nrow(outInfo), " records.")
    # message("Downloaded  ", round(i/nrow(genesURL)*100,2),"%; totally ", length(na.omit(outInfo$gencodeId)), " records fetched!")
    # rm(url1, url1Get, url1GetText, url1GetText2Json, url1GetText2Json2DT)
  }
  outInfo <- merge(outInfo, tissueSiteDetailGTEx, by="tissueSiteDetailId")
  outInfo <- merge(outInfo[,.(gencodeId, geneSymbol, median, tissueSiteDetail)], geneInfo[,.(gencodeId, geneSymbol, entrezGeneId, geneType, chromosome, start, end, strand, tss)], by=c("gencodeId", "geneSymbol"))
  message("Unit of expression: ",unique(tmp$unit))
  return(outInfo)
}

#' @title Retrive SNP pairwise LD from locuscompare database.
#' @description
#'  SNP pairwise lD are calculated based on 1000 Genomes Project Phase 3 version 5.
#'  For storage-efficiency, the output will only include SNPs with r2 > 0.2 with the input SNP.
#' @param chr (string) Chromosome name. e.g. '22'. Notice that the name should not contain 'chr'.
#' @param snp (string) SNP rsID.
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @import RMySQL
#' @import DBI
#' @export
#' @examples
#' \donttest{
#'  ld <- retrieveLD('6', 'rs9349379', 'AFR')
#'  }
retrieveLD = function(chr,snp,population){
  # conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
  conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare", "locuscomparer" ,"12345678","locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com")
  on.exit(RMySQL::dbDisconnect(conn))

  res1 = DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_A, SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_A = '%s';",
      chr,
      population,
      snp
    )
  )

  res2 = DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_B as SNP_A, SNP_A as SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_B = '%s';",
      chr,
      population,
      snp
    )
  )

  res = rbind(res1,res2)
  return(res)
}

#' @title Retrive SNP pairwise LD from LDlink database
#'
#' @param targetSnp target SNP, support dbSNP IP.
#' @param population Supported population is consistent with the LDlink, which can be listed using function LDlinkR::list_pop()
#' @param windowSize Window around the highlighted snp for querying linkage disequilibrium information. Default:500000
#' @param method The same as fetchContent function, can be chosen from "download", "curl", "GetWithHeader", or "GET".
#' @param genomeVersion "grch38"(default) or "grch37".
#' @param max_count To prevent download failure due to network fluctuations, max number of connection attempts.
#' @param token Ldlink token, default: "9246d2db7917"
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' \donttest{
#'   snpLD <- retrieveLD_LDproxy("rs3", windowSize=5000)
#' }
retrieveLD_LDproxy <- function(targetSnp="", population="EUR" , windowSize=500000, method="download", genomeVersion="grch38", max_count=3, token="9246d2db7917"){
  # targetSnp="rs3"
  # population="EUR"
  # windowSize=500000
  # genomeVersion="grch38"
  # max_count=3
  # token="9246d2db7917"

  # get LD information:
  s_count<-1
  max_count<- 3

  while( !exists("snpLDtmp") && s_count<=max_count ){
    message("=== Geting LD info for SNP: ",targetSnp,"; trying ",s_count,"/",max_count,".")
    url1 <- paste0("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=",targetSnp,
                   "&pop=",population,
                   "&r2_d=","r2",
                   "&window=",as.character(as.integer(windowSize)),
                   "&genome_build=",genomeVersion,
                   "&token=", token)

    # url1 <- "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=rs3&pop=MXL&r2_d=r2&window=100000&genome_build=grch38&token=9246d2db7917"
    try( snpLDtmp <- fetchContent(url1, method="download", isJson=FALSE) )
    if( exists("snpLDtmp") && ncol(snpLDtmp)<=1 ){
      rm(snpLDtmp)
    }else{
      message("=== Done!")
    }
    s_count <- s_count+1
  }
  if(!exists("snpLDtmp") ){
    return( NULL )
  }else{
    return(snpLDtmp)
  }
}

#' @title fetch matrix of pairwise linkage disequilibrium statistics.
#' @description This fucntion is rebuilt from LDlink package.
#' @param snps list of between 2 - 1,000 variants, using an rsID or chromosome coordinate (e.g. "chr7:24966446")
#' @param pop a 1000 Genomes Project population, (e.g. YRI or CEU), multiple allowed, default = "CEU". All supported can be listed with LDlinkR::list_pop()
#' @param r2d r2d, either "r2" for LD R2 or "d" for LD D', default = "r2"
#' @param token LDlink provided user token, default = NULL, register for token at https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param file Optional character string naming a path and file for saving results. If file = FALSE, no file will be generated, default = FALSE.
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \donttest{
#'   ldMat <- retrieveLD_LDmatrix(c("rs12202891", "rs10807323","rs9381401", "rs2026458", "rs150503442"),
#'                       pop="CEU", token="9246d2db7917")
#' }
retrieveLD_LDmatrix <- function(snps, pop = "CEU", r2d = "r2", token = NULL,  file = FALSE) {
  LD_config <- list(ldmatrix_url = "https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix",
                    avail_pop = c("YRI", "LWK", "GWD",
                                  "MSL", "ESN", "ASW", "ACB",
                                  "MXL", "PUR", "CLM", "PEL",
                                  "CHB", "JPT", "CHS", "CDX",
                                  "KHV", "CEU", "TSI", "FIN",
                                  "GBR", "IBS", "GIH", "PJL",
                                  "BEB", "STU", "ITU", "ALL",
                                  "AFR", "AMR", "EAS", "EUR",
                                  "SAS"), avail_ld = c("r2", "d"))
  url <- LD_config[["ldmatrix_url"]]
  avail_pop <- LD_config[["avail_pop"]]
  avail_ld <- LD_config[["avail_ld"]]
  file <- as.character(file)
  rsid_pattern <- "^rs\\d{1,}"
  chr_coord_pattern <- "(^chr)(\\d{1,2}|X|x|Y|y):(\\d{1,9})$"
  if (!(length(snps) > 1) & (length(snps) <= 1000)) {
    stop("Input is between 2 to 1000 variants.")
  }
  for (i in 1:length(snps)) {
    if (!((grepl(rsid_pattern, snps[i], ignore.case = TRUE)) |
          (grepl(chr_coord_pattern, snps[i], ignore.case = TRUE)))) {
      stop(paste("Invalid query format for variant: ",
                 snps[i], ".", sep = ""))
    }
  }
  if (!(all(pop %in% avail_pop))) {
    stop("Not a valid population code.")
  }
  if (!(r2d %in% avail_ld)) {
    stop("Not a valid r2d.  Enter 'r2' or 'd'.")
  }
  if (is.null(token)) {
    stop("Enter valid access token. Please register using the LDlink API Access tab: https://ldlink.nci.nih.gov/?tab=apiaccess")
  }
  if (!(is.character(file) | file == FALSE)) {
    stop("Invalid input for file option.")
  }
  snps_to_upload <- paste(unlist(snps), collapse = "\n")
  pop_to_upload <- paste(unlist(pop), collapse = "+")
  jsonbody <- list(snps = snps_to_upload, pop = pop_to_upload,
                   r2_d = r2d)
  url_str <- paste(url, "?", "&token=", token,
                   sep = "")
  # message(url_str)
  if (httr::http_error(url)) {
    message("The LDlink server is down or not accessible. Please try again later.")
    return(NULL)
  }
  else {
    message("\nLDlink server is working...\n")
  }
  raw_out <- httr::POST(url = url_str, body = jsonbody, encode = "json")
  httr::stop_for_status(raw_out)
  data_out <- read.delim(textConnection(httr::content(raw_out,
                                                      "text", encoding = "UTF-8")), header = T,
                         sep = "\t")
  if (grepl("error", data_out[2, 1])) {
    stop(data_out[2, 1])
  }
  if (file == FALSE) {
    return(data_out)
  }
  else if (is.character(file)) {
    print(data_out)
    write.table(data_out, file = file, quote = F, row.names = F,
                sep = "\t")
    cat(paste("\nFile saved to ", file, ".",
              sep = ""))
    return(data_out)
  }
}

