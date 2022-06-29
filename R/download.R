#' @title Download normalized gene expression at the sample level in a specified tissue.
#' @param genes (character string or a character vector) gene symbols or gencode ids (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @param toSummarizedExperiment a logical value indicating whether to return a data.frame or a summarizedExperiment object. Default: TRUE, return a toSummarizedExperiment object.
#' @param recordPerChunk (integer) number of records fetched per request (default: 80).
#' @param pathologyNotesCategories a logical value indicating whether to return pathologyNotes. Default: FALSE, the pathologyNotes is ignored.
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import stats
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom   IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @return return a SummarizedExperiment or a data.table object harbing gene expression profiles and samples' information.
#' @export
#' @examples
#' \donttest{
#' # Download gene expression with a genecode ID:
#' expProfiles <- xQTLdownload_exp("ENSG00000210195.2", tissueSiteDetail="Liver")
#'
#' # extract expression profile from SummarizedExperiment object:
#' expDT <- SummarizedExperiment::assay(expProfiles)
#'
#' # extract samples' detail from SummarizedExperiment object:
#' sampleDT <- SummarizedExperiment::colData(expProfiles)
#'
#' # Download gene expression profiles of multiple genes:
#' expProfiles <- xQTLdownload_exp(c("tp53","naDK","SDF4"),
#'                                 tissueSiteDetail="Artery - Coronary",
#'                                 pathologyNotesCategories=TRUE,
#'                                 toSummarizedExperiment=FALSE)
#'
#' # Download with versioned and unversioned gencode Id.
#' expProfiles <- xQTLdownload_exp(c("ENSG00000141510.16","ENSG00000008130.15","ENSG00000078808"),
#'                                tissueSiteDetail="Artery - Coronary",
#'                                toSummarizedExperiment=FALSE)
#' }
xQTLdownload_exp <- function(genes="", geneType="auto", tissueSiteDetail="Liver", datasetId="gtex_v8", toSummarizedExperiment=TRUE, recordPerChunk=80, pathologyNotesCategories=FALSE  ){
  gencodeId <- chromosome <- cutF <- genesUpper <- geneSymbol <- entrezGeneId <- tss <- description <- NULL
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
    if( all(unlist(lapply(genes, function(g){ stringr::str_detect(g, "^ENSG") }))) ){
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
  message("== Check the gene name :")
  geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion)
  # Only keep genes with non-na gencode ID
  geneInfo <- geneInfo[!is.na(gencodeId)]
  # for gene of the same name but with different gencodeID, like SHOX, ENSG00000185960.13 the name is in X, ENSG00000185960.13_PAR_Y is the name in Y. retain the gencode In x.
  geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & stringr::str_detect(chromosome, "Y"))]
  # for gene of the same name but with different gencodeID, like BTBD8 with entrezGeneId is NA, remove.
  geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & is.na(entrezGeneId))]
  # retain source is Source:NCBI
  geneInfo[(genes %in% geneInfo[duplicated(geneSymbol), ]$genes)]

  # duplicates, only warning first duplicate:
  # test: genes = c("LYNX1", "TP53")
  # dupGeneSymbol <- unique(geneInfo[,.(geneSymbol, gencodeId)][,.(.SD[duplicated(geneSymbol)])]$geneSymbol[1])
  # if(!is.na(dupGeneSymbol)){
  #   message(" ")
  #   message("  == Gene [",dupGeneSymbol,"] has multiple gencode IDs: [",paste0(geneInfo[geneSymbol == dupGeneSymbol]$gencodeId, collapse = ", "),"]")
  #   stop("== Please remove duplicated genes, or take the unique gencode IDs as the input.")
  # }
  dupGencodeId <- unique(geneInfo[,.(geneSymbol, gencodeId)][,.(.SD[duplicated(gencodeId)])]$gencodeId[1])
  if(!is.na(dupGencodeId)){
    stop("duplicated gencodeId detected: ", paste0(dupGencodeId, collapse = "; "))
  }

  #
  if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
    stop("gene information is null.")
  }
  message("== Done.")
  # geneInfo <-  xQTLquery_gene(c("tp53","naDK","SDF4"), "geneSymbol", "v26", "GRCh38/hg38")

  ############ get sample info:
  message("== Get the samples' detail:")
  sampleInfo <- xQTLquery_sampleByTissue(tissueSiteDetail=tissueSiteDetail, dataType="RNASEQ", datasetId=datasetId, pathologyNotesCategories=pathologyNotesCategories )
  message("== Done.")
  if( !exists("sampleInfo") ||is.null(sampleInfo) ){
    stop("Failed to fetch sample information.")
  }

  ############### 分批下载：
  message("== Get the gene expression:")
  geneInfo <- cbind(geneInfo, data.table::data.table(cutF = cut(1:nrow(geneInfo),breaks=seq(0,nrow(geneInfo)+cutNum,cutNum) ) ))
  data.table::setDT(geneInfo)
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
    tmp <- merge(geneInfo[cutF==genesURL[i,]$cutF,.(genes, gencodeId, geneSymbol, datasetId)], tmp[,-c("geneSymbol", "datasetId")], by="gencodeId",all.x=TRUE, sort = FALSE)
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
  rownames(expDT) <- outInfo$genes
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


#' @title Download significant eQTL associations of a specified tissue or across all tissues.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param genes (character string or a character vector) gene symbols or gencode ids (versioned or unversioned are both supported).
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' # Download eQTL info for a variant:
#' xQTLdownload_eqtlSig("rs201327123")
#' xQTLdownload_eqtlSig("chr1_14677_G_A_b38")
#' xQTLdownload_eqtlSig("11_66328719_T_C_b37", datasetId="gtex_v7")
#' xQTLdownload_eqtlSig("11_66328719_T_C_b37", datasetId="gtex_v7",
#'                     tissueSiteDetail="Skin - Sun Exposed (Lower leg)")
#'
#' # Download eQTL association according to all tissues with genome location:
#' varInfo <-  xQTLquery_varPos(chrom="chr1", pos=c(1102708),"gtex_v8")
#' xQTLdownload_eqtlSig(variantName=varInfo$snpId)
#'
#' # Download eQTL info for gene:
#' xQTLdownload_eqtlSig(genes="ATAD3B")
#' xQTLdownload_eqtlSig(genes=c("TP53", "SLC35E2B"), tissueSiteDetail= "Brain - Cerebellum")
#' xQTLdownload_eqtlSig(genes="ENSG00000141510.16", datasetId="gtex_v8")
#' xQTLdownload_eqtlSig(genes="ENSG00000141510.11", datasetId="gtex_v7",
#'                      tissueSiteDetail="Thyroid" )
#'
#' # Download eQTL info for a variant-gene pair:
#' xQTLdownload_eqtlSig(variantName="rs1641513", genes="TP53", datasetId="gtex_v8")
#' xQTLdownload_eqtlSig(variantName="rs1641513", genes="TP53", datasetId="gtex_v7")
#' xQTLdownload_eqtlSig(variantName="chr1_1667948_A_G_b38",
#'                      genes="SLC35E2B", tissueSiteDetail="Kidney - Cortex")
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

    message("== Check the variant name:")
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

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=genes, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    geneInfo <- geneInfo[!is.na(gencodeId)]
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


#' @title Download significant or unsignificant eQTL associations of a tissue or across all tissues
#' @description
#'  can be quried with a gene/variant-gene pair.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported). Can not be null.
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param recordPerChunk (integer) number of records fetched per request (default: 100).
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' # Download eQTL info with a gene symbol:
#' eqtlInfo <- xQTLdownload_eqtl(gene="TP53")
#'
#' # Use unversioned gencode ID in GTEx V8:
#' eqtl_v8 <- xQTLdownload_eqtl(gene="ENSG00000141510", datasetId="gtex_v8")
#'
#' # In a specific tissue:
#' xQTLdownload_eqtl(gene="ENSG00000141510.11", geneType="gencodeId",
#'                   datasetId="gtex_v7", tissueSiteDetail="Thyroid" )
#'
#' # Download eQTL info with a variant-gene pair:
#' xQTLdownload_eqtl(variantName="rs1641513",gene="TP53", datasetId="gtex_v8")
#' xQTLdownload_eqtl(variantName="chr1_1667948_A_G_b38", gene="SLC35E2B",
#'                   tissueSiteDetail="Kidney - Cortex")
#' }
xQTLdownload_eqtl <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail="", datasetId="gtex_v8", recordPerChunk=100){
  pos <- variantId <- snpId <- gencodeId <- geneSymbol <- pValue <- nes <- se <- mValue <- NULL
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

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    geneInfo <- geneInfo[!is.na(gencodeId)]
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

    message("== Check the variant name:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input, or leave \"variantName\" as null.")
    }else{
      message("== Done.")
    }
  }

  # url1 <- "https://gtexportal.org/rest/v1/association/metasoft?gencodeId=ENSG00000141510.11&datasetId=gtex_v7"
  message("== Querying eQTL associations from API server:")
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


#' @title Download summary statistics of eQTL of a specified gene, variant, tissue or study.
#' @description
#'  source of all eQTL associations is EBI eQTL category.
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueLabel (character) all supported tissues can be listed using "ebi_study_tissues"
#' @param recordPerChunk (integer) number of records fetched per request (default: 1000).
#' @param study (character) name of studies can be listed using "ebi_study_tissues". If the study is null, use all studies (Default).
#' @param withB37VariantId a logical value indicating whether to return the genome location(GTEx v7) of variants. Default: FALSE.
#' @import data.table
#' @import stringr
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' # Download all associations of MLH1-rs13315355 pair in all tissues from all studies:
#' eqtlAsso <- xQTLdownload_eqtlAllAsso(gene="MLH1", variantName = "rs13315355", study="")
#'
#' # Download associations of gene ATP11B in CD4+ T cell from all supported studies(time-consuming):
#' geneAsso <- xQTLdownload_eqtlAllAsso(gene="MMP7",tissueLabel = "CD4+ T cell", study="")
#'
#' # Download associations of gene ATP11B in Muscle - Skeletal from GTEx_V8:
#' geneAsso <- xQTLdownload_eqtlAllAsso("ATP11B", tissueLabel="Muscle - Skeletal")
#'
#' # Download all associations of SNP rs11568818 in all tissues from all supported studies.
#' varAsso <- xQTLdownload_eqtlAllAsso(variantName="rs11568818", study="")
#
#' # Download associations of SNP rs11568818 in Muscle - Skeletal from GTEx_V8:
#' varAsso <- xQTLdownload_eqtlAllAsso(variantName="chr11_102530930_T_C_b38",
#'                                     tissueLabel="Muscle - Skeletal")
#' }
xQTLdownload_eqtlAllAsso <- function(gene="", geneType="auto", variantName="", variantType="auto", tissueLabel="", study="gtex_v8", recordPerChunk=1000, withB37VariantId=FALSE){
  . <- geneInfoV19 <- pos <- ref<- alt<- tissue<-NULL
  chrom <- tissue_label <- study_accession <- variantId <- variant <- gencodeId <- genes<- entrezGeneId <- chromosome<- geneSymbol<- b37VariantId <- snpId <- NULL
  # gene="CYP2W1"
  # geneType="geneSymbol"
  # tissueSiteDetail="Lung"

  ebi_ST <-copy(ebi_study_tissues)

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

  # check tissue:
  if( length(tissueLabel)==1 && tissueLabel!="" ){
    if( toupper(tissueLabel) %in% toupper(unique(ebi_ST$tissue_label)) ){
      message("== Tissue label [", tissueLabel, "] detected...")
      tissueLabel <- unique(ebi_ST$tissue_label)[ toupper(unique(ebi_ST$tissue_label)) == toupper(tissueLabel) ]
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== tissueLabel [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
    }
  }

  # check study-tissue:
  if( length(study) ==1 && length(tissueLabel)==1 && study!="" && tissueLabel!=""){
    if(nrow( ebi_ST[study_accession == study & tissue_label==tissueLabel])==1){
      message("== Study [",study,"] -- Tissue label [",tissueLabel,"] corrected mapped..")
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== Study [",study,"] -- Tissue label [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
    }
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
  # modify variant name of variantID:
  if(length(variantName)==1 && variantName!="" && variantType == "variantId"){
    variantName <- stringr::str_remove(variantName, "_b38")
  }

  # Excessive empty variables:
  if(gene=="" && variantName=="" && tissueLabel==""){
    warning("== All associations from Study [", study,"] will be returned, this will take several days...")
    warning("== Please make sure you have enough memory...")
    message("== Specified gene, variant or tissueLabel is recommended.")
  }
  if(gene=="" && variantName=="" && study==""){
    warning("== All associations in Tissue [", tissueLabel,"] will be returned, this will take several days...")
    warning("== Please make sure you have enough memory...")
    message("== Specified gene, variant or study is recommended.")
  }
  if(gene=="" && variantName=="" && tissueLabel=="" && study==""){
    stop("== All required fileds are empty, please specify gene, variant study or tissueLabel")
  }

  # if gene is not null, check geneName and add gene unversioned ensemble name.
  if(gene !=""){
    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v26")
    # 对于只需要一个基因的：
    geneInfo <- geneInfo[!is.na(gencodeId)]
    geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & stringr::str_detect(chromosome, "Y"))]
    geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & is.na(entrezGeneId))]
    # geneInfoV19 <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v19")
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or set gene with gencodeId.")
    }else{
      geneInfo$gencodeIdUnv <-stringr::str_split(geneInfo$gencodeId, stringr::fixed("."))[[1]][1]
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if( variantName!="" && variantType!= "auto"){
    # construct url:
    # aa<-"https://www.ebi.ac.uk/eqtl/api/associations/rs2302765?study=GTEx_V8&gene_id=ENSG00000238917&tissue=UBER_0001211"
    url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/associations/",variantName,"?links=False",
                   ifelse(gene !="", paste0("&gene_id=",geneInfo$gencodeIdUnv),""),
                   ifelse(tissueLabel!="", paste0("&tissue=",ebi_ST[tissue_label==tissueLabel]$tissue[1]),""),
                   ifelse(study!="", paste0("&study=",study),"")
                   )
  }else if(gene!="" && geneType!="auto"){
    # aa <-“https://www.ebi.ac.uk/eqtl/api/genes/ENSG00000282458/associations?variant_id=chr19_80901_G_T&tissue=UBER_0001954&study=GTEx_V8”
    url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/genes/",geneInfo$gencodeIdUnv,"/associations?links=False",
                   ifelse(tissueLabel!="", paste0("&tissue=",ebi_ST[tissue_label==tissueLabel]$tissue[1]),""),
                   ifelse(study!="", paste0("&study=",study),"")
                   # ifelse(variantName!="",paste0("&variant_id=",variantName),"")
                   )
  }else if( tissueLabel!=""  ){
    # aa <- "https://www.ebi.ac.uk/eqtl/api/tissues/CL_0000235/studies/Alasoo_2018/associations"
    url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/tissues/",ebi_ST[tissue_label==tissueLabel]$tissue[1],
                   ifelse(study!="", paste0("/studies/",study),""),
                   "/associations?links=False"
                   )
  }else if( study!="" ){
    # aa <- "https://www.ebi.ac.uk/eqtl/api/studies/GTEx_V8/associations"
    url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/studies/", study, "/associations?links=False")
  }

  # for brain tissues with duplicated tissue id, "Brain - Cerebellar Hemisphere" and "Brain - Cerebellum"
  if( tissueLabel == "Brain - Cerebellar Hemisphere" ){
    qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellar_Hemisphere")
  }else if(tissueLabel == "Brain - Cerebellum"){
    qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellum")
  }else{
    qtl_groupStr <- ""
  }
  url1 <- paste0(url1, qtl_groupStr)
  # message(url1)

  message("== Start fetching associations...")
  # construct url:

  # # check network:
  # bestFetchMethod <- apiEbi_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: EBI API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # gtexAsoo <- fetchContentEbi(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2],  termSize=1000)
  gtexAsoo <- fetchContentEbi(url1, method = "fromJSON", downloadMethod = "auto",  termSize=recordPerChunk)
  if(is.null(gtexAsoo)){
    return(NULL)
  }
  gtexAsooList <- do.call(c, gtexAsoo)

  if(length(gtexAsooList)==0){
    message("No association found!")
    return(NULL)
  }
  gtexAsooList <- lapply(gtexAsooList, function(x){ data.table(variantId= x$variant, snpId=x$rsid,type=x$type,maf=x$maf,beta=x$beta,
                                                               chrom=x$chromosome, pos=x$position,ref=x$ref, alt=x$alt,study_id=x$study_id,
                                                               se=x$se, median_tpm=x$median_tpm, pValue=x$pvalue, totalAlleles=x$an, allelCounts=x$ac, imputationR2=x$r2,
                                                               tissue_label=x$tissue_label,tissue=x$tissue,qtl_group=x$qtl_group, condition=x$condition,
                                                               molecular_trait_id = x$molecular_trait_id, gene_id= x$gene_id
                                                               ) })
  gtexAsooDT <- data.table::rbindlist(gtexAsooList, fill=TRUE)
  gtexAsooDT[,variantId:=.(paste0( "chr",stringr::str_remove_all(chrom, "chr"),"_",pos,"_",ref,"_",alt ))]
  if( !("tissue_label" %in% names(gtexAsooDT)) ){
    if(tissueLabel==""){
      ebi_ST_tmp <- ebi_ST[tissue_label!="Brain - Cerebellum", .(tissue, tissue_label,study_id=study_accession)]
      }else{ ebi_ST_tmp <-ebi_ST[tissue_label==tissueLabel,.(tissue, tissue_label,study_id=study_accession)]  }
    gtexAsooDT <- merge(gtexAsooDT, ebi_ST_tmp, by = c("tissue","study_id"))
  }
  gtexAsooDT <- cbind(gtexAsooDT[,-c("study_id", "tissue", "tissue_label", "qtl_group")], gtexAsooDT[,c("study_id", "tissue", "tissue_label", "qtl_group")])
  if(gene!="" && geneType!="auto"){
    gtexAsooDT$geneSymbol <- geneInfo$geneSymbol
    gtexAsooDT$gencodeId_GTEX_v8 <- geneInfo$gencodeId
    # gtexAsooDT$gencodeId_GTEX_v7 <- ifelse( exists("geneInfoV19") && nrow(geneInfoV19)>0, geneInfoV19$gencodeId, "")
  }else{
    geneInfo <- xQTLquery_gene(unique(gtexAsooDT$gene_id))
    gtexAsooDT <- merge(gtexAsooDT[,], geneInfo[,.(gene_id=genes, geneSymbol, gencodeId_GTEX_v8=gencodeId)], by = "gene_id", all.x = TRUE)[,-c("gene_id")]
  }
  if(withB37VariantId){
    # add dbSNP id and  hg19 cordinate:
    gtexAsooDTb37 <- xQTLquery_varPos(chrom = paste0("chr",unique(gtexAsooDT$chrom)), pos = unique(gtexAsooDT$pos), datasetId = "gtex_v8", recordPerChunk = 250)
    gtexAsooDTb37$variantId <- unlist(lapply(gtexAsooDTb37$variantId, function(x){ splitInfo=stringr::str_split(x, stringr::fixed("_"))[[1]]; paste0(splitInfo[-5], collapse="_") }))
    gtexAsooDT <- merge(gtexAsooDT, gtexAsooDTb37[,.(variantId, b37VariantId)], by=c("variantId"), all.x=TRUE )
    # gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
    # gtexAsooDT <- cbind(gtexAsooDT[,.(snpId, variantId, b37VariantId)], gtexAsooDT[,-c("snpId", "variantId", "b37VariantId", "chrom", "pos")])
  }
  gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
  return(gtexAsooDT)
}



#' @title Download summary statistics of eQTL with genome position.
#'
#' @param chrom (character) name of chromesome, including chr1-chr22, chrX.
#' @param pos_lower (integer) lower base pair location threshold, expressed as an integer
#' @param pos_upper (integer) upper base pair location threshold, expressed as an integer
#' @param p_lower (numeric) lower p-value threshold, can be expressed as a float or using mantissa and exponent annotation (0.001 or 1e-3 or 1E-3)
#' @param p_upper (numeric) upper p-value threshold, can be expressed as a float or using mantissa and exponent annotation (0.001 or 1e-3 or 1E-3)
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueLabel (character) all supported tissues can be listed using "ebi_study_tissues".
#' @param study (character) name of studies can be listed using "ebi_study_tissues".
#' @param recordPerChunk (integer) number of records fetched per request (default: 1000).
#' @param withB37VariantId a logical value indicating whether to return the genome location(GTEx v7) of variants. Default: FALSE.
#'
#' @import data.table
#' @import stringr
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#' eqtlAssos <- xQTLdownload_eqtlAllAssoPos(chrom = "chr11",
#'                                         pos_lower=101398614, pos_upper = 101402313,
#'                                         tissueLabel="Brain - Cerebellar Hemisphere",
#'                                         p_upper=1e-1)
#' }
xQTLdownload_eqtlAllAssoPos <- function(chrom="", pos_lower=numeric(0), pos_upper=numeric(0), p_lower=0, p_upper=1.1,  gene="", geneType="auto", tissueLabel="", study="gtex_v8", recordPerChunk=1000, withB37VariantId=FALSE){
  .<-NULL
  study_accession <- tissue_label <- gencodeId <- genes <- geneSymbol <- chromosome <- entrezGeneId <- variantId <- pos <- ref <- alt <- tissue <- b37VariantId <- NULL
  ebi_ST <- data.table::copy(ebi_study_tissues)

  # check chrome:
  if(length(chrom)==1 && chrom!="" && !is.na(chrom) ){
    chrom <- toupper(stringr::str_remove(chrom, "chr"))
  }else if( !(chrom %in% c(1:22, "X")) ){
    stop("Chromosome must be choosen from one of: ", paste0(paste0("chr",c(1:22, "X")), collapse = ", "))
  }else{
    message("== chromosome: ", chrom)
  }

  # check position and pvalue:
  if( length(pos_lower)==1 && length(pos_upper) ==1 && length(p_upper)==1 && length(p_lower)==1 && is.wholenumber(pos_lower) && is.wholenumber(pos_upper) && is.numeric(p_lower) && is.numeric(p_upper)){
    message("== Genome range: ",chrom,":", as.character(as.integer(pos_lower)),"-",as.character(as.integer(pos_upper)))
  }else{
    stop("pos_lower, pos_upper must be the whole numbers, and p_lower, p_upper must be the numeric value.")
  }
  if(pos_upper < pos_lower){stop("pos_upper must greater than pos_lower")}
  if( p_upper < p_lower ){ stop("p_upper must greater than p_lower")}


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

  # check tissue:
  if( length(tissueLabel)==1 && tissueLabel!="" ){
    if( toupper(tissueLabel) %in% toupper(unique(ebi_ST$tissue_label)) ){
      message("== Tissue label [", tissueLabel, "] detected...")
      tissueLabel <- unique(ebi_ST$tissue_label)[ toupper(unique(ebi_ST$tissue_label)) == toupper(tissueLabel) ]
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== tissueLabel [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
    }
  }

  # check study-tissue:
  if( length(study) ==1 && length(tissueLabel)==1 && study!="" && tissueLabel!=""){
    if(nrow( ebi_ST[study_accession == study & tissue_label==tissueLabel])==1){
      message("== Study [",study,"] -- Tissue label [",tissueLabel,"] corrected mapped..")
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== Study [",study,"] -- Tissue label [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
    }
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

  # if gene is not null, check geneName and add gene unversioned ensemble name.
  if(gene !=""){
    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v26")
    # 对于只需要一个基因的：
    geneInfo <- geneInfo[!is.na(gencodeId)]
    geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & stringr::str_detect(chromosome, "Y"))]
    geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & is.na(entrezGeneId))]
    # geneInfoV19 <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v19")
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or set gene with gencodeId.")
    }else{
      geneInfo$gencodeIdUnv <-stringr::str_split(geneInfo$gencodeId, stringr::fixed("."))[[1]][1]
      message("== Done.")
    }
  }

  message("== Start fetching associations...", format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))
  # url1 <- "https://www.ebi.ac.uk/eqtl/api/chromosomes/11/associations?study=GTEx_V8&gene_id=ENSG00000137673&bp_lower=101798614&bp_upper=103462313&p_upper=1e-1"
  url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/chromosomes/", chrom, "/associations?links=False",
                 "&bp_lower=",as.character(as.integer(pos_lower)), "&bp_upper=",as.character(as.integer(pos_upper)), "&p_upper=",p_upper, "&p_lower=",p_lower,
                 ifelse(gene !="", paste0("&gene_id=",geneInfo$gencodeIdUnv),""),
                 ifelse(tissueLabel!="", paste0("&tissue=",ebi_ST[tissue_label==tissueLabel]$tissue[1]),""),
                 ifelse(study!="", paste0("&study=",study),"")
                 )

  # for brain tissues with duplicated tissue id, "Brain - Cerebellar Hemisphere" and "Brain - Cerebellum"
  if( tissueLabel == "Brain - Cerebellar Hemisphere" ){
    qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellar_Hemisphere")
  }else if(tissueLabel == "Brain - Cerebellum"){
    qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellum")
  }else{
    qtl_groupStr <- ""
  }
  url1 <- paste0(url1, qtl_groupStr)

  # message(url1)
  gtexAsoo <- fetchContentEbi(url1, method = "fromJSON", downloadMethod = "auto",  termSize=recordPerChunk)
  gtexAsooList <- do.call(c, gtexAsoo)
  if(length(gtexAsooList)==0){
    message("No association found!")
    return(NULL)
  }
  gtexAsooList <- lapply(gtexAsooList, function(x){ data.table(variantId= x$variant, snpId=x$rsid,type=x$type,maf=x$maf,beta=x$beta,
                                                               chrom=x$chromosome, pos=x$position,ref=x$ref, alt=x$alt,study_id=x$study_id,
                                                               se=x$se, median_tpm=x$median_tpm, pValue=x$pvalue, totalAlleles=x$an, allelCounts=x$ac, imputationR2=x$r2,
                                                               tissue_label=x$tissue_label,tissue=x$tissue,qtl_group=x$qtl_group, condition=x$condition,
                                                               molecular_trait_id = x$molecular_trait_id, gene_id= x$gene_id
  ) })
  gtexAsooDT <- data.table::rbindlist(gtexAsooList, fill=TRUE)
  gtexAsooDT[,variantId:=.(paste0( "chr",stringr::str_remove_all(chrom, "chr"),"_",pos,"_",ref,"_",alt ))]
  if( !("tissue_label" %in% names(gtexAsooDT)) ){
    if(tissueLabel==""){
      ebi_ST_tmp <- ebi_ST[tissue_label!="Brain - Cerebellum", .(tissue, tissue_label,study_id=study_accession)]
    }else{ ebi_ST_tmp <-ebi_ST[tissue_label==tissueLabel,.(tissue, tissue_label,study_id=study_accession)]  }
    gtexAsooDT <- merge(gtexAsooDT, ebi_ST_tmp, by = c("tissue","study_id"))
  }
  gtexAsooDT <- cbind(gtexAsooDT[,-c("study_id", "tissue", "tissue_label", "qtl_group")], gtexAsooDT[,c("study_id", "tissue", "tissue_label", "qtl_group")])

  # add gene info.
  if(gene!="" && geneType!="auto"){
    gtexAsooDT$geneSymbol <- geneInfo$geneSymbol
    gtexAsooDT$gencodeId_GTEX_v8 <- geneInfo$gencodeId
    # gtexAsooDT$gencodeId_GTEX_v7 <- ifelse( exists("geneInfoV19") && nrow(geneInfoV19)>0, geneInfoV19$gencodeId, "")
  }else{
    geneInfo <- xQTLquery_gene(unique(gtexAsooDT$gene_id))
    gtexAsooDT <- merge(gtexAsooDT[,], geneInfo[,.(gene_id=genes, geneSymbol, gencodeId_GTEX_v8=gencodeId)], by = "gene_id", all.x = TRUE)[,-c("gene_id")]
  }
  if(withB37VariantId){
    # add dbSNP id and  hg19 cordinate:
    gtexAsooDTb37 <- xQTLquery_varPos(chrom = paste0("chr",unique(gtexAsooDT$chrom)), pos = unique(gtexAsooDT$pos), datasetId = "gtex_v8", recordPerChunk = 250)
    gtexAsooDTb37$variantId <- unlist(lapply(gtexAsooDTb37$variantId, function(x){ splitInfo=stringr::str_split(x, stringr::fixed("_"))[[1]]; paste0(splitInfo[-5], collapse="_") }))
    gtexAsooDT <- merge(gtexAsooDT, gtexAsooDTb37[,.(variantId, b37VariantId)], by=c("variantId"), all.x=TRUE )
    # gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
    # gtexAsooDT <- cbind(gtexAsooDT[,.(snpId, variantId, b37VariantId)], gtexAsooDT[,-c("snpId", "variantId", "b37VariantId", "chrom", "pos")])
  }
  gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
  return(gtexAsooDT)
}


#' @title Download significant sQTL associations of a tissue or across all tissues
#' @description
#'  Only GTEx v8 is supported.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param genes (character string or a character vector) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' # Download sQTL details with rsid:
#' xQTLdownload_sqtlSig(variantName="rs201327123")
#' xQTLdownload_sqtlSig(variantName="chr1_14677_G_A_b38", tissueSiteDetail="Whole Blood")
#'
#' # Download sQTL details with gene symbol and gencode ID:
#' xQTLdownload_sqtlSig(genes="ENSG00000141510.16", tissueSiteDetail="Lung" )
#' xQTLdownload_sqtlSig(genes=c("ATAD3B", "MLH1"))
#'
#' # Download sQTL details with the variant-genes pair:
#' xQTLdownload_sqtlSig(variantName="rs201327123", genes=c("WASH7P","RP11-206L10.2"))
#' xQTLdownload_sqtlSig(variantName="chr17_7465085_A_G_b38", genes="TP53", tissueSiteDetail="Lung")
#' }
xQTLdownload_sqtlSig <- function(variantName="", genes="", variantType="auto", geneType="auto", tissueSiteDetail="" ){
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

    message("== Check the variant name:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if( all(genes !="")){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=genes, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    geneInfo <- geneInfo[!is.na(gencodeId)]
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or leave \"gene\" as null")
    }else{
      message("== Done.")
    }
  }
  # both are null: error!
  if(variantName=="" && all(genes == "") ){
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
                 ifelse(all(genes==""),"",paste0("&gencodeId=",paste0(geneInfo$gencodeId, collapse = ","))),
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
  if(!exists("tmp")||nrow(tmp)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & genes!="","-",""),ifelse(genes=="","",paste0(" gene [", genes,"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail)), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  tmp <- merge(tmp, tissueSiteDetailGTEx, by = "tissueSiteDetailId")
  outInfo <- tmp[,.(variantId, snpId, gencodeId, geneSymbol, phenotypeId, tissueSiteDetail, pValue, nes, datasetId)]
  message("=================================")
  message("Totally ", nrow(outInfo), " associatons were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & genes!=""," -",""),ifelse(all(genes==""),"",paste0(" gene: [",paste0(genes, collapse = ", "),"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail))," in ",datasetId,"."  )

  return(outInfo)
}

#' @title Download normalized expression of gene for a eQTL pair.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import stats
#' @return A data.table object.
#' @export
#' @examples
#' \donttest{
#' # Download exp with variant-gene pair in different tissues:
#' xQTLdownload_eqtlExp(variantName="rs1641513",gene="TP53", tissueSiteDetail="Liver")
#'
#' # Download expression using variant ID and gencode ID.
#' xQTLdownload_eqtlExp(variantName="chr1_14677_G_A_b38",gene="ENSG00000228463.9",
#'                      tissueSiteDetail="Stomach")
#' }
xQTLdownload_eqtlExp <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail="", datasetId="gtex_v8"){
  gencodeGenetype <- chromosome <-gencodeId <-NULL
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

    message("== Check the variant name:")
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

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    geneInfo <- geneInfo[!is.na(gencodeId)]
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


#' @title Download normalized expression of intron for a sQTL pair.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param phenotypeId A character string. Format like: "chr1:497299:498399:clu_54863:ENSG00000239906.1"
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import stats
#' @return A data.table object.
#' @export
#'
#' @examples
#' # Download sQTL expression in different tissues:
#' xQTLdownload_sqtlExp(variantName="rs1450891501",
#'                      phenotypeId="chr1:497299:498399:clu_54863:ENSG00000239906.1",
#'                      tissueSiteDetail="Lung")
#'
#' # Dowload sQTL expression using variant ID.
#' xQTLdownload_sqtlExp(variantName="chr1_1259424_T_C_b38",
#'                      phenotypeId=" chr1:1487914:1489204:clu_52051:ENSG00000160072.19",
#'                      tissueSiteDetail="Adipose - Subcutaneous")
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

    message("== Check the variant name:")
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



#' @title Download details of eGenes (eQTL Genes) for a specified gene or a tissue.
#' @description
#'  eGenes are genes that have at least one significant cis-eQTL acting upon them. Results can be filtered by tissue.
#' @param gene  (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param recordPerChunk (integer) number of records fetched per request (default: 2000).
#' @import data.table
#' @import stringr
#' @import utils
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' eGeneInfo <- xQTLdownload_egene("TP53")
#' eGeneInfo <- xQTLdownload_egene(tissueSiteDetail="Prostate", recordPerChunk=2000)
#' eGeneInfo <- xQTLdownload_egene("DDX11", datasetId="gtex_v7",tissueSiteDetail="Artery - Tibial")
#' }
xQTLdownload_egene <- function(gene = "", geneType="auto", datasetId = "gtex_v8", tissueSiteDetail="", recordPerChunk=2000){
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

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType, gencodeVersion=gencodeVersion)
    geneInfo <- geneInfo[!is.na(gencodeId)]
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


#' @title Download details of sGenes (sQTL Genes) for a specified gene or a tissue.
#' @description
#'  sGenes are genes that have at least one significant sQTL acting upon them. Results may be filtered by tissue.
#' @param gene  (character) gene symbol or gencode id (versioned or unversioned are both supported). Can be null.
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param recordPerChunk (integer) number of records fetched per request (default: 2000).
#' @import data.table
#' @import stringr
#' @import utils
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' sGeneInfo <- xQTLdownload_sgene(tissueSiteDetail="Liver")
#' sGeneInfo <- xQTLdownload_sgene(gene="DDX11", tissueSiteDetail="Liver" )
#' }
xQTLdownload_sgene <- function(gene = "", geneType="auto", datasetId = "gtex_v8", tissueSiteDetail="", recordPerChunk=2000){
  phenotypeId <- nPhenotypes <- .<-NULL
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

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType, gencodeVersion=gencodeVersion)
    geneInfo <- geneInfo[!is.na(gencodeId)]
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


#' @title Download median expression of all samples for specified genes across tissues.
#' @param genes (character string or a character vector) gene symbols or gencode ids (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param recordPerChunk (integer) number of records fetched per request (default: 150).
#' @import data.table
#' @import utils
#' @import stringr
#' @return A data.table object.
#' @export
#'
#' @examples
#' geneMedExp <- xQTLdownload_geneMedExp(genes="LYNX1")
#' geneMedExp <- xQTLdownload_geneMedExp(genes=c("TP53", "IRF5"))
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
  if(all(genes!="") || length(genes)>0){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Check the gene name:")
    suppressMessages( geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, recordPerChunk = recordPerChunk))
    geneInfo <- geneInfo[!is.na(gencodeId)]
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
#' @return A data.frame object.
#' @export
#' @examples
#' ld <- retrieveLD('6', 'rs9349379', 'AFR')
retrieveLD = function(chr, snp, population){
  # conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
  conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare", "locuscomparer" ,"12345678","locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com")
  on.exit(RMySQL::dbDisconnect(conn))

  chr <- str_remove_all(chr,"chr")
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






