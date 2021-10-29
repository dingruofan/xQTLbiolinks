
#' @title Fetch normalized gene expression.
#' @description
#'  This function fetch queried genes' expression profiles in a tissue at the sample level.
#' @param genes Following gene types are supported:
#' \itemize{
#'   \item \strong{Gene symbol}.
#'
#'   A character string or a character vector (case ignored). like: "tp53","naDK","SDF4".
#'   \item \strong{Gencode/ensemble id} (versioned or unversioned).
#' }
#'
#' @param geneType A character string. Types of queried genes. Options: "geneSymbol" (default) and "gencodeId";
#' @param tissueSiteDetail
#'  Tissue must be chosen from the following tissue names:
#' \tabular{rrrrr}{
#'   \strong{tissue name} \tab \strong{GTEx V8} \tab \strong{GTEx V7} \cr
#'    Adipose - Subcutaneous \tab √ \tab √\cr
#'    Adipose - Visceral (Omentum) \tab √ \tab √\cr
#'    Adrenal Gland \tab √ \tab √\cr
#'    Artery - Aorta \tab √ \tab √\cr
#'    Artery - Coronary \tab √ \tab √\cr
#'    Artery - Tibial \tab √ \tab √\cr
#'    Bladder \tab √ \tab √\cr
#'    Brain - Amygdala \tab √ \tab √\cr
#'    Brain - Anterior cingulate cortex (BA24) \tab √ \tab √\cr
#'    Brain - Caudate (basal ganglia) \tab √ \tab √\cr
#'    Brain - Cerebellar Hemisphere \tab √ \tab √\cr
#'    Brain - Cerebellum \tab √ \tab √\cr
#'    Brain - Cortex \tab √ \tab √\cr
#'    Brain - Frontal Cortex (BA9) \tab √ \tab √\cr
#'    Brain - Hippocampus \tab √ \tab √\cr
#'    Brain - Hypothalamus \tab √ \tab √\cr
#'    Brain - Nucleus accumbens (basal ganglia) \tab √ \tab √\cr
#'    Brain - Putamen (basal ganglia) \tab √ \tab √\cr
#'    Brain - Spinal cord (cervical c-1) \tab √ \tab √\cr
#'    Brain - Substantia nigra \tab √ \tab √\cr
#'    Breast - Mammary Tissue \tab √ \tab √\cr
#'    Cells - Cultured fibroblasts \tab √ \tab x\cr
#'    Cells - EBV-transformed lymphocytes \tab √ \tab √\cr
#'    Cells - Transformed fibroblasts \tab x \tab √\cr
#'    Cervix - Ectocervix \tab √ \tab √\cr
#'    Cervix - Endocervix \tab √ \tab √\cr
#'    Colon - Sigmoid \tab √ \tab √\cr
#'    Colon - Transverse \tab √ \tab √\cr
#'    Esophagus - Gastroesophageal Junction \tab √ \tab √\cr
#'    Esophagus - Mucosa \tab √ \tab √\cr
#'    Esophagus - Muscularis \tab √ \tab √\cr
#'    Fallopian Tube \tab √ \tab √\cr
#'    Heart - Atrial Appendage \tab √ \tab √\cr
#'    Heart - Left Ventricle \tab √ \tab √\cr
#'    Kidney - Cortex \tab √ \tab √\cr
#'    Kidney - Medulla \tab √ \tab x\cr
#'    Liver \tab √ \tab √\cr
#'    Lung \tab √ \tab √\cr
#'    Minor Salivary Gland \tab √ \tab √\cr
#'    Muscle - Skeletal \tab √ \tab √\cr
#'    Nerve - Tibial \tab √ \tab √\cr
#'    Ovary \tab √ \tab √\cr
#'    Pancreas \tab √ \tab √\cr
#'    Pituitary \tab √ \tab √\cr
#'    Prostate \tab √ \tab √\cr
#'    Skin - Not Sun Exposed (Suprapubic) \tab √ \tab √\cr
#'    Skin - Sun Exposed (Lower leg) \tab √ \tab √\cr
#'    Small Intestine - Terminal Ileum \tab √ \tab √\cr
#'    Spleen \tab √ \tab √\cr
#'    Stomach \tab √ \tab √\cr
#'    Testis \tab √ \tab √\cr
#'    Thyroid \tab √ \tab √\cr
#'    Uterus \tab √ \tab √\cr
#'    Vagina \tab √ \tab √\cr
#'    Whole Blood \tab √ \tab √\cr
#' }
#' @param datasetId A character string. Options: "gtex_v8" (default), "gtex_v7".
#' @param toSummarizedExperiment whether to return a data.frame or a summarizedExperiment object. Default: TRUE, return a toSummarizedExperiment object.
#' @param recordPerChunk A integer value (1-2000). number of records fetched per request (default: 150).
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
#'   expProfiles <- GTExdownload_exp(c("ENSG00000210195.2"),
#'                                  "gencodeId", "Liver", "gtex_v8")
#'   # Download gene expression profiles with multiple genes:
#'   expProfiles <- GTExdownload_exp(c("tp53","naDK","SDF4"),
#'                                  "geneSymbol", "Artery - Coronary", "gtex_v8",
#'                                  toSummarizedExperiment=TRUE)
#'   expProfiles <- GTExdownload_exp(c("tp53","naDK","SDF4"),
#'                                  "geneSymbol", "Artery - Coronary", "gtex_v7")
#'
#'
#'   # Get proteing-coding genes' expression:
#'   proT <- GTExquery_gene (genes="protein coding", geneType="geneCategory", "v26" )
#'   proTexp <- GTExdownload_exp(proT$geneSymbol[1:100], geneType = "geneSymbol","Lung","gtex_v8")
#'   }
GTExdownload_exp <- function(genes="", geneType="geneSymbol", tissueSiteDetail="Liver", datasetId="gtex_v8", toSummarizedExperiment=TRUE, recordPerChunk=150  ){
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

  ############ convert genes. parameter check is unnecessary for this, because GTExquery_gene check it internally.
  message("== Querying gene info from API server:")
  geneInfo <- GTExquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, recordPerChunk=recordPerChunk)
  if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
    stop("gene information is null.")
  }
  message("== Done.")
  # geneInfo <-  GTExquery_gene(c("tp53","naDK","SDF4"), "geneSymbol", "v26", "GRCh38/hg38")

  ############ get sample info:
  message("== Fetching sample information from API server:")
  sampleInfo <- GTExquery_sample(tissueSiteDetail=tissueSiteDetail, dataType="RNASEQ", datasetId=datasetId, recordPerChunk=recordPerChunk )
  message("== Done.")
  if( is.null(sampleInfo)||!exists("sampleInfo") ){
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
  pingOut <- apiAdmin_ping()
  if( !(!is.null(pingOut) && pingOut==200) ){
    return(data.table::data.table())
  }
  for(i in 1:nrow(genesURL)){
    # construct url:
    url1 <- paste0("https://gtexportal.org/rest/v1/expression/geneExpression?",
                   "datasetId=", datasetId,"&",
                   "gencodeId=", genesURL[i,]$genesURL, "&",
                   "tissueSiteDetailId=", tissueSiteDetailId,"&",
                   "&format=json"
    )
    url1 <- utils::URLencode(url1)
    url1Get <- curl::curl_fetch_memory(url1)
    if(url1Get$status_code!=200){
      stop("Http status code: ", url1Get$status_code)
    }
    url1GetText <- rawToChar(url1Get$content)
    url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
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
    rm(url1, url1Get, url1GetText, url1GetText2Json, tmp)
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

#' @title Download significant eQTL data.
#' @description
#'  Fetch significant eQTL associations with a variant or a gene or a variant-gene pair in a tissue or across all tissues.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "snpId" or "variantId". Default: "snpId".
#' @param geneType A character string. "geneSymbol"(default) or "gencodeId".
#' @param tissueSiteDetail A character string.
#'  Tissue must be chosen from the following tissue names:
#' \tabular{rrrrr}{
#'   \strong{tissue name} \tab \strong{GTEx V8} \tab \strong{GTEx V7} \cr
#'    Adipose - Subcutaneous \tab √ \tab √\cr
#'    Adipose - Visceral (Omentum) \tab √ \tab √\cr
#'    Adrenal Gland \tab √ \tab √\cr
#'    Artery - Aorta \tab √ \tab √\cr
#'    Artery - Coronary \tab √ \tab √\cr
#'    Artery - Tibial \tab √ \tab √\cr
#'    Bladder \tab √ \tab √\cr
#'    Brain - Amygdala \tab √ \tab √\cr
#'    Brain - Anterior cingulate cortex (BA24) \tab √ \tab √\cr
#'    Brain - Caudate (basal ganglia) \tab √ \tab √\cr
#'    Brain - Cerebellar Hemisphere \tab √ \tab √\cr
#'    Brain - Cerebellum \tab √ \tab √\cr
#'    Brain - Cortex \tab √ \tab √\cr
#'    Brain - Frontal Cortex (BA9) \tab √ \tab √\cr
#'    Brain - Hippocampus \tab √ \tab √\cr
#'    Brain - Hypothalamus \tab √ \tab √\cr
#'    Brain - Nucleus accumbens (basal ganglia) \tab √ \tab √\cr
#'    Brain - Putamen (basal ganglia) \tab √ \tab √\cr
#'    Brain - Spinal cord (cervical c-1) \tab √ \tab √\cr
#'    Brain - Substantia nigra \tab √ \tab √\cr
#'    Breast - Mammary Tissue \tab √ \tab √\cr
#'    Cells - Cultured fibroblasts \tab √ \tab x\cr
#'    Cells - EBV-transformed lymphocytes \tab √ \tab √\cr
#'    Cells - Transformed fibroblasts \tab x \tab √\cr
#'    Cervix - Ectocervix \tab √ \tab √\cr
#'    Cervix - Endocervix \tab √ \tab √\cr
#'    Colon - Sigmoid \tab √ \tab √\cr
#'    Colon - Transverse \tab √ \tab √\cr
#'    Esophagus - Gastroesophageal Junction \tab √ \tab √\cr
#'    Esophagus - Mucosa \tab √ \tab √\cr
#'    Esophagus - Muscularis \tab √ \tab √\cr
#'    Fallopian Tube \tab √ \tab √\cr
#'    Heart - Atrial Appendage \tab √ \tab √\cr
#'    Heart - Left Ventricle \tab √ \tab √\cr
#'    Kidney - Cortex \tab √ \tab √\cr
#'    Kidney - Medulla \tab √ \tab x\cr
#'    Liver \tab √ \tab √\cr
#'    Lung \tab √ \tab √\cr
#'    Minor Salivary Gland \tab √ \tab √\cr
#'    Muscle - Skeletal \tab √ \tab √\cr
#'    Nerve - Tibial \tab √ \tab √\cr
#'    Ovary \tab √ \tab √\cr
#'    Pancreas \tab √ \tab √\cr
#'    Pituitary \tab √ \tab √\cr
#'    Prostate \tab √ \tab √\cr
#'    Skin - Not Sun Exposed (Suprapubic) \tab √ \tab √\cr
#'    Skin - Sun Exposed (Lower leg) \tab √ \tab √\cr
#'    Small Intestine - Terminal Ileum \tab √ \tab √\cr
#'    Spleen \tab √ \tab √\cr
#'    Stomach \tab √ \tab √\cr
#'    Testis \tab √ \tab √\cr
#'    Thyroid \tab √ \tab √\cr
#'    Uterus \tab √ \tab √\cr
#'    Vagina \tab √ \tab √\cr
#'    Whole Blood \tab √ \tab √\cr
#' }
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
#'  GTExdownload_eqtlSig(variantName="rs201327123", variantType="snpId")
#'  GTExdownload_eqtlSig(variantName="chr1_14677_G_A_b38", variantType="variantId")
#'  GTExdownload_eqtlSig(variantName="11_66328719_T_C_b37", variantType="variantId",
#'                    datasetId="gtex_v7")
#'  GTExdownload_eqtlSig(variantName="11_66328719_T_C_b37", variantType="variantId",
#'                    datasetId="gtex_v7", tissueSiteDetail="Skin - Sun Exposed (Lower leg)")
#'
#'  # Download eQTL association according to all tissues with genome location:
#'  varInfo <-  GTExquery_varPos(chrom="chr1", pos=c(1102708),"gtex_v8")
#'  GTExdownload_eqtlSig(variantName=varInfo$snpId, variantType="snpId")
#'
#'  # Download eQTL info for a gene:
#'  GTExdownload_eqtlSig(gene="ATAD3B")
#'  GTExdownload_eqtlSig(gene="TP53",datasetId="gtex_v7")
#'  GTExdownload_eqtlSig(gene="ENSG00000141510.16", geneType="gencodeId", datasetId="gtex_v8")
#'  GTExdownload_eqtlSig(gene="ENSG00000141510.11", geneType="gencodeId",
#'                    datasetId="gtex_v7",tissueSiteDetail="Thyroid" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  GTExdownload_eqtlSig(variantName="rs1641513",gene="TP53", datasetId="gtex_v8")
#'  GTExdownload_eqtlSig(variantName="rs1641513",gene="TP53", datasetId="gtex_v7")
#'  GTExdownload_eqtlSig(variantName="chr1_1667948_A_G_b38", variantType="variantId",
#'                    gene="SLC35E2B", geneType="geneSymbol", tissueSiteDetail="Kidney - Cortex")
#' }
GTExdownload_eqtlSig <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8"){
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
  pingOut <- apiAdmin_ping()
  if( !(!is.null(pingOut) && pingOut==200) ){
    return(data.table::data.table())
  }
  ##################### fetch varInfo
  if(variantName!=""){
    message("== Querying variant info from API server:")
    varInfo <- GTExquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input, or leave \"variantName\" as null.")
    }else{
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if(gene !=""){
    message("== Querying gene info from API server:")
    geneInfo <- GTExquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
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
  message("== Querying significant eQTL associations from API server:")
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/singleTissueEqtl?format=json",
                 "&datasetId=",datasetId,
                 ifelse(variantName=="","",paste0("&variantId=",varInfo$variantId)),
                 ifelse(gene=="","",paste0("&gencodeId=",geneInfo$gencodeId)),
                 ifelse(tissueSiteDetail=="","",paste0("&tissueSiteDetailId=",tissueSiteDetailId))
                 )
  url1 <- utils::URLencode(url1)
  url1Get <- curl::curl_fetch_memory(url1)
  if(url1Get$status_code!=200){
    stop("Http status code: ", url1Get$status_code)
  }
  url1GetText <- rawToChar(url1Get$content)
  url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
  tmp <- data.table::as.data.table(url1GetText2Json$singleTissueEqtl)
  if(nrow(tmp)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & gene!="","-",""),ifelse(gene=="","",paste0(" gene [", gene,"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail)), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  tmp <- merge(tmp, tissueSiteDetailGTEx, by = "tissueSiteDetailId")
  outInfo <- tmp[,.(variantId, snpId, gencodeId, geneSymbol, tissueSiteDetail, pValue, nes,datasetId)]
  message("=================================")
  message("Totally ", nrow(outInfo), " associatons were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & gene!=""," -",""),ifelse(gene=="","",paste0(" gene [", gene,"]")),ifelse(tissueSiteDetail=="",paste0(" in ",length(unique(tmp$tissueSiteDetail)),ifelse(length(unique(tmp$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail))," in ",datasetId,"."  )

  return(outInfo)
}


#' @title Download significant or unsignificant eQTL data.
#' @description
#'  Fetch significant or unsignificant eQTL associations with a gene or a variant-gene pair in a tissue or across all tissues.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned). Can not be null.
#' @param variantType A character string. "snpId" or "variantId". Default: "snpId".
#' @param geneType A character string. "geneSymbol"(default) or "gencodeId".
#' @param tissueSiteDetail A character string.
#'  Tissue must be chosen from the following tissue names:
#' \tabular{rrrrr}{
#'   \strong{tissue name} \tab \strong{GTEx V8} \tab \strong{GTEx V7} \cr
#'    Adipose - Subcutaneous \tab √ \tab √\cr
#'    Adipose - Visceral (Omentum) \tab √ \tab √\cr
#'    Adrenal Gland \tab √ \tab √\cr
#'    Artery - Aorta \tab √ \tab √\cr
#'    Artery - Coronary \tab √ \tab √\cr
#'    Artery - Tibial \tab √ \tab √\cr
#'    Bladder \tab √ \tab √\cr
#'    Brain - Amygdala \tab √ \tab √\cr
#'    Brain - Anterior cingulate cortex (BA24) \tab √ \tab √\cr
#'    Brain - Caudate (basal ganglia) \tab √ \tab √\cr
#'    Brain - Cerebellar Hemisphere \tab √ \tab √\cr
#'    Brain - Cerebellum \tab √ \tab √\cr
#'    Brain - Cortex \tab √ \tab √\cr
#'    Brain - Frontal Cortex (BA9) \tab √ \tab √\cr
#'    Brain - Hippocampus \tab √ \tab √\cr
#'    Brain - Hypothalamus \tab √ \tab √\cr
#'    Brain - Nucleus accumbens (basal ganglia) \tab √ \tab √\cr
#'    Brain - Putamen (basal ganglia) \tab √ \tab √\cr
#'    Brain - Spinal cord (cervical c-1) \tab √ \tab √\cr
#'    Brain - Substantia nigra \tab √ \tab √\cr
#'    Breast - Mammary Tissue \tab √ \tab √\cr
#'    Cells - Cultured fibroblasts \tab √ \tab x\cr
#'    Cells - EBV-transformed lymphocytes \tab √ \tab √\cr
#'    Cells - Transformed fibroblasts \tab x \tab √\cr
#'    Cervix - Ectocervix \tab √ \tab √\cr
#'    Cervix - Endocervix \tab √ \tab √\cr
#'    Colon - Sigmoid \tab √ \tab √\cr
#'    Colon - Transverse \tab √ \tab √\cr
#'    Esophagus - Gastroesophageal Junction \tab √ \tab √\cr
#'    Esophagus - Mucosa \tab √ \tab √\cr
#'    Esophagus - Muscularis \tab √ \tab √\cr
#'    Fallopian Tube \tab √ \tab √\cr
#'    Heart - Atrial Appendage \tab √ \tab √\cr
#'    Heart - Left Ventricle \tab √ \tab √\cr
#'    Kidney - Cortex \tab √ \tab √\cr
#'    Kidney - Medulla \tab √ \tab x\cr
#'    Liver \tab √ \tab √\cr
#'    Lung \tab √ \tab √\cr
#'    Minor Salivary Gland \tab √ \tab √\cr
#'    Muscle - Skeletal \tab √ \tab √\cr
#'    Nerve - Tibial \tab √ \tab √\cr
#'    Ovary \tab √ \tab √\cr
#'    Pancreas \tab √ \tab √\cr
#'    Pituitary \tab √ \tab √\cr
#'    Prostate \tab √ \tab √\cr
#'    Skin - Not Sun Exposed (Suprapubic) \tab √ \tab √\cr
#'    Skin - Sun Exposed (Lower leg) \tab √ \tab √\cr
#'    Small Intestine - Terminal Ileum \tab √ \tab √\cr
#'    Spleen \tab √ \tab √\cr
#'    Stomach \tab √ \tab √\cr
#'    Testis \tab √ \tab √\cr
#'    Thyroid \tab √ \tab √\cr
#'    Uterus \tab √ \tab √\cr
#'    Vagina \tab √ \tab √\cr
#'    Whole Blood \tab √ \tab √\cr
#' }
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
#'  GTExdownload_eqtlAll(gene="TP53")
#'  GTExdownload_eqtlAll(gene="ATAD3B", datasetId="gtex_v7")
#'
#'  # Unversioned gencode ID in GTEx V7:
#'  eqtl_v7 <- GTExdownload_eqtlAll(gene="ENSG00000141510", geneType="gencodeId",
#'                                  datasetId="gtex_v7")
#'  # Unversioned gencode ID in GTEx V8:
#'  eqtl_v8 <- GTExdownload_eqtlAll(gene="ENSG00000141510", geneType="gencodeId",
#'                                  datasetId="gtex_v8")
#'  # No intersection of eQTLs for gene ENSG00000141510!
#'  intersect(eqtl_v7$snpId, eqtl_v8$snpId)
#'
#'  # In a specific tissue:
#'  GTExdownload_eqtlAll(gene="ENSG00000141510.11", geneType="gencodeId",
#'                       datasetId="gtex_v7", tissueSiteDetail="Thyroid" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  GTExdownload_eqtlAll(variantName="rs1641513",gene="TP53", datasetId="gtex_v8")
#'  GTExdownload_eqtlAll(variantName="rs1641513",gene="TP53",
#'                       datasetId="gtex_v7")
#'  GTExdownload_eqtlSig(variantName="chr1_1667948_A_G_b38", variantType="variantId",
#'                       gene="SLC35E2B", geneType="geneSymbol",
#'                       tissueSiteDetail="Kidney - Cortex")
#'  GTExdownload_eqtlAll(variantName="chr1_1667948_A_G_b38", variantType="variantId",
#'                       gene="SLC35E2B", geneType="geneSymbol",
#'                       tissueSiteDetail="Kidney - Cortex")
#' }
GTExdownload_eqtlAll <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8"){
  pos <- variantId <- snpId <- gencodeId <- geneSymbol <- pValue <- nes <- se <- NULL
  .<-NULL
  cutNum<-100
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
  pingOut <- apiAdmin_ping()
  if( !(!is.null(pingOut) && pingOut==200) ){
    return(data.table::data.table())
  }
  ##################### fetch geneInfo:
  if(gene !=""){
    message("== Querying gene info from API server:")
    geneInfo <- GTExquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
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
    message("== Querying variant info from API server:")
    varInfo <- GTExquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
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
  url1Get <- curl::curl_fetch_memory(url1)
  if(url1Get$status_code!=200){
    stop("Http status code: ", url1Get$status_code)
  }
  url1GetText <- rawToChar(url1Get$content)
  url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
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
    tmp_i <- cbind( tmp_info, tmp[,names(tmp)[i]])
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
    pingOut <- apiAdmin_ping()
    if( !(!is.null(pingOut) && pingOut==200) ){
      return(data.table::data.table())
    }
    var_tmp <- data.table::data.table(variantId =unique(outInfo$variantId))

    ##### query with GTExquery_varId: START
    # var_tmp$snpId <- unlist(lapply(1:nrow(var_tmp), function(x){message("   Got varints: ",x,"/",nrow(var_tmp));GTExquery_varId(var_tmp$variantId[x], variantType = "variantId", datasetId = datasetId)$snpId; }))
    ##### query with GTExquery_varId: END

    ##### query with GTExquery_varPos: START
    var_tmp <- cbind(var_tmp, data.table::rbindlist(lapply(unique(outInfo$variantId), function(x){ splitOut = stringr::str_split(x,stringr::fixed("_"))[[1]];data.table::data.table(chrom=splitOut[1], pos=splitOut[2]) })))
    # query pos in batch per 100 terms.
    var_tmpCut <- cbind(var_tmp, data.table::data.table( ID=1:nrow(var_tmp), cutF = as.character(cut(1:nrow(var_tmp),breaks=seq(0,nrow(var_tmp)+cutNum,cutNum) )) ))
    var_tmpCut <- var_tmpCut[,.(posLis=list(pos)),by=c("chrom","cutF")]
    # out:
    var_tmpGot <- data.table::data.table()
    for( ii in 1:nrow(var_tmpCut)){
      message("   Got varints: ",nrow(var_tmpGot)+length(unlist(var_tmpCut[ii,]$posLis)),"/",nrow(var_tmp))
      var_tmpGot <- rbind(var_tmpGot, GTExquery_varPos(var_tmpCut[ii,]$chrom, pos=as.integer(unlist(var_tmpCut[ii,]$posLis)), datasetId = datasetId))
    }
    var_tmp <- merge(var_tmp[,.(variantId)], var_tmpGot, by="variantId", all.x=TRUE)
    ##### query with GTExquery_varPos: END

    outInfo <- merge(outInfo, var_tmp, by="variantId")
    message("== Done.")
  }else{
    outInfo$snpId <- ""
  }
  outInfo <- outInfo[,.(variantId, snpId, gencodeId, geneSymbol, tissueSiteDetail, pValue, nes, se, datasetId)]

  message("=================================")
  message("In total of ", nrow(outInfo) ," eQTL associations were found for ", ifelse(gene=="","",paste0("Gene [", gene,"]")), ifelse(variantName=="","",paste0(" - Variant [",variantName,"]")), ifelse(tissueSiteDetail=="", paste0(" in ",length(unique(outInfo$tissueSiteDetail)),ifelse(length(unique(outInfo$tissueSiteDetail))==1," tissue", " tissues")), paste0(" in ", tissueSiteDetail))," in ",datasetId,"."  )
  return(outInfo)
}


#' @title Download significant sQTL data.
#' @description
#'  Fetch significant sQTL associations with a variant or a gene or a variant-gene pair in a tissue or across all tissues.
#'  Only GTEx v8 is supported.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "snpId" or "variantId". Default: "snpId".
#' @param geneType A character string. "geneSymbol"(default) or "gencodeId".
#' @param tissueSiteDetail A character string.
#'  Tissue must be chosen from the following tissue names:
#' \tabular{rrrrr}{
#'   \strong{tissue name} \tab \strong{GTEx V8} \cr
#'    Adipose - Subcutaneous \tab √ \cr
#'    Adipose - Visceral (Omentum) \tab √ \cr
#'    Adrenal Gland \tab √ \cr
#'    Artery - Aorta \tab √ \cr
#'    Artery - Coronary \tab √ \cr
#'    Artery - Tibial \tab √ \cr
#'    Bladder \tab √ \cr
#'    Brain - Amygdala \tab √ \cr
#'    Brain - Anterior cingulate cortex (BA24) \tab √ \cr
#'    Brain - Caudate (basal ganglia) \tab √ \cr
#'    Brain - Cerebellar Hemisphere \tab √ \cr
#'    Brain - Cerebellum \tab √ \cr
#'    Brain - Cortex \tab √ \cr
#'    Brain - Frontal Cortex (BA9) \tab √ \cr
#'    Brain - Hippocampus \tab √ \cr
#'    Brain - Hypothalamus \tab √ \cr
#'    Brain - Nucleus accumbens (basal ganglia) \tab √ \cr
#'    Brain - Putamen (basal ganglia) \tab √ \cr
#'    Brain - Spinal cord (cervical c-1) \tab √ \cr
#'    Brain - Substantia nigra \tab √ \cr
#'    Breast - Mammary Tissue \tab √ \cr
#'    Cells - Cultured fibroblasts \tab √ \cr
#'    Cells - EBV-transformed lymphocytes \tab √ \cr
#'    Cells - Transformed fibroblasts \tab x \cr
#'    Cervix - Ectocervix \tab √ \cr
#'    Cervix - Endocervix \tab √ \cr
#'    Colon - Sigmoid \tab √ \cr
#'    Colon - Transverse \tab √ \cr
#'    Esophagus - Gastroesophageal Junction \tab √ \cr
#'    Esophagus - Mucosa \tab √ \cr
#'    Esophagus - Muscularis \tab √ \cr
#'    Fallopian Tube \tab √ \cr
#'    Heart - Atrial Appendage \tab √ \cr
#'    Heart - Left Ventricle \tab √ \cr
#'    Kidney - Cortex \tab √ \cr
#'    Kidney - Medulla \tab √ \cr
#'    Liver \tab √ \cr
#'    Lung \tab √ \cr
#'    Minor Salivary Gland \tab √ \cr
#'    Muscle - Skeletal \tab √ \cr
#'    Nerve - Tibial \tab √ \cr
#'    Ovary \tab √ \cr
#'    Pancreas \tab √ \cr
#'    Pituitary \tab √ \cr
#'    Prostate \tab √ \cr
#'    Skin - Not Sun Exposed (Suprapubic) \tab √ \cr
#'    Skin - Sun Exposed (Lower leg) \tab √ \cr
#'    Small Intestine - Terminal Ileum \tab √ \cr
#'    Spleen \tab √ \cr
#'    Stomach \tab √ \cr
#'    Testis \tab √ \cr
#'    Thyroid \tab √ \cr
#'    Uterus \tab √ \cr
#'    Vagina \tab √ \cr
#'    Whole Blood \tab √ \cr
#' }
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
#'  GTExdownload_sqtlSig(variantName="rs201327123", variantType="snpId")
#'  GTExdownload_sqtlSig(variantName="chr1_14677_G_A_b38", variantType="variantId")
#'  GTExdownload_sqtlSig(variantName="chr1_14677_G_A_b38", variantType="variantId",
#'                       tissueSiteDetail="Whole Blood")
#'
#'  # Download sQTL association according to all tissues with genome location:
#'  varInfo <-  GTExquery_varPos(chrom="chr1", pos=c(1102708),"gtex_v8")
#'  GTExdownload_sqtlSig(variantName=varInfo$snpId, variantType="snpId")
#'
#'  # Download sQTL info for a gene:
#'  GTExdownload_sqtlSig(gene="ATAD3B")
#'  GTExdownload_sqtlSig(gene="ENSG00000141510.16", geneType="gencodeId")
#'  GTExdownload_sqtlSig(gene="ENSG00000141510.16", geneType="gencodeId",
#'                       tissueSiteDetail="Lung" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  GTExdownload_sqtlSig(variantName="rs546057177", gene="TP53")
#'  GTExdownload_sqtlSig(variantName="chr17_7465085_A_G_b38", variantType="variantId",
#'                       gene="TP53", geneType="geneSymbol",
#'                       tissueSiteDetail="Lung")
#' }
GTExdownload_sqtlSig <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="" ){
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
  # check network:
  pingOut <- apiAdmin_ping()
  if( !(!is.null(pingOut) && pingOut==200) ){
    return(data.table::data.table())
  }
  ##################### fetch varInfo
  if( variantName!=""){
    message("== Querying variant info from API server:")
    varInfo <- GTExquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if(gene !=""){
    message("== Querying gene info from API server:")
    geneInfo <- GTExquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
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
  url1 <- utils::URLencode(url1)
  url1Get <- curl::curl_fetch_memory(url1)
  if(url1Get$status_code!=200){
    stop("Http status code: ", url1Get$status_code)
  }
  url1GetText <- rawToChar(url1Get$content)
  url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
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

#' @title Download expression data for eQTL.
#' @description
#'  This function fetch normalized gene expression data for a eQTL pair.
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "snpId" or "variantId". Default: "snpId".
#' @param geneType A character string. "geneSymbol"(default) or "gencodeId".
#' @param tissueSiteDetail A character string.
#'  Tissue must be chosen from the following tissue names:
#' \tabular{rrrrr}{
#'   \strong{tissue name} \tab \strong{GTEx V8} \tab \strong{GTEx V7} \cr
#'    Adipose - Subcutaneous \tab √ \tab √\cr
#'    Adipose - Visceral (Omentum) \tab √ \tab √\cr
#'    Adrenal Gland \tab √ \tab √\cr
#'    Artery - Aorta \tab √ \tab √\cr
#'    Artery - Coronary \tab √ \tab √\cr
#'    Artery - Tibial \tab √ \tab √\cr
#'    Bladder \tab √ \tab √\cr
#'    Brain - Amygdala \tab √ \tab √\cr
#'    Brain - Anterior cingulate cortex (BA24) \tab √ \tab √\cr
#'    Brain - Caudate (basal ganglia) \tab √ \tab √\cr
#'    Brain - Cerebellar Hemisphere \tab √ \tab √\cr
#'    Brain - Cerebellum \tab √ \tab √\cr
#'    Brain - Cortex \tab √ \tab √\cr
#'    Brain - Frontal Cortex (BA9) \tab √ \tab √\cr
#'    Brain - Hippocampus \tab √ \tab √\cr
#'    Brain - Hypothalamus \tab √ \tab √\cr
#'    Brain - Nucleus accumbens (basal ganglia) \tab √ \tab √\cr
#'    Brain - Putamen (basal ganglia) \tab √ \tab √\cr
#'    Brain - Spinal cord (cervical c-1) \tab √ \tab √\cr
#'    Brain - Substantia nigra \tab √ \tab √\cr
#'    Breast - Mammary Tissue \tab √ \tab √\cr
#'    Cells - Cultured fibroblasts \tab √ \tab x\cr
#'    Cells - EBV-transformed lymphocytes \tab √ \tab √\cr
#'    Cervix - Ectocervix \tab √ \tab √\cr
#'    Cervix - Endocervix \tab √ \tab √\cr
#'    Colon - Sigmoid \tab √ \tab √\cr
#'    Colon - Transverse \tab √ \tab √\cr
#'    Esophagus - Gastroesophageal Junction \tab √ \tab √\cr
#'    Esophagus - Mucosa \tab √ \tab √\cr
#'    Esophagus - Muscularis \tab √ \tab √\cr
#'    Fallopian Tube \tab √ \tab √\cr
#'    Heart - Atrial Appendage \tab √ \tab √\cr
#'    Heart - Left Ventricle \tab √ \tab √\cr
#'    Kidney - Cortex \tab √ \tab √\cr
#'    Kidney - Medulla \tab √ \tab x\cr
#'    Liver \tab √ \tab √\cr
#'    Lung \tab √ \tab √\cr
#'    Minor Salivary Gland \tab √ \tab √\cr
#'    Muscle - Skeletal \tab √ \tab √\cr
#'    Nerve - Tibial \tab √ \tab √\cr
#'    Ovary \tab √ \tab √\cr
#'    Pancreas \tab √ \tab √\cr
#'    Pituitary \tab √ \tab √\cr
#'    Prostate \tab √ \tab √\cr
#'    Skin - Not Sun Exposed (Suprapubic) \tab √ \tab √\cr
#'    Skin - Sun Exposed (Lower leg) \tab √ \tab √\cr
#'    Small Intestine - Terminal Ileum \tab √ \tab √\cr
#'    Spleen \tab √ \tab √\cr
#'    Stomach \tab √ \tab √\cr
#'    Testis \tab √ \tab √\cr
#'    Thyroid \tab √ \tab √\cr
#'    Uterus \tab √ \tab √\cr
#'    Vagina \tab √ \tab √\cr
#'    Whole Blood \tab √ \tab √\cr
#' }
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import stats
#' @return A data.table
#' @export
#' @examples
#' \donttest{
#'  # Download exp in different tissues:
#'  GTExdownload_eqtlExp(variantName="rs1641513",gene="TP53", tissueSiteDetail="Liver")
#'  GTExdownload_eqtlExp(variantName="rs1641513",gene="ATAD3B",
#'                       tissueSiteDetail="Lung", datasetId="gtex_v8")
#'
#'  # Download exp in gtex v7:
#'  GTExdownload_eqtlExp(variantName="rs140894808",gene="ATAD3B",
#'                       tissueSiteDetail="Adipose - Visceral (Omentum)",
#'                       datasetId="gtex_v7")
#'
#'  # Dowload exp using variant ID and gencode ID.
#'  GTExdownload_eqtlExp(variantName="chr1_14677_G_A_b38",gene="ENSG00000228463.9",
#'                       variantType="variantId", geneType="gencodeId",
#'                      tissueSiteDetail="Stomach")
#' }
GTExdownload_eqtlExp <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8"){
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

  # check network:
  pingOut <- apiAdmin_ping()
  if( !(!is.null(pingOut) && pingOut==200) ){
    return(data.table::data.table())
  }
  ##################### fetch varInfo
  if(variantName!=""){
    message("== Querying variant info from API server:")
    varInfo <- GTExquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
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
    message("== Querying gene info from API server:")
    geneInfo <- GTExquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    if(nrow(geneInfo)==0|| is.null(geneInfo) || !exists("geneInfo")){
      stop("Invalid gene name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"gene\" can not be null!")
  }

  # construct url
  message("== Querying expression from API server:")
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/dyneqtl?",
                 "gencodeId=",geneInfo$gencodeId,"&",
                 "variantId=",varInfo$variantId,"&",
                 "tissueSiteDetailId=",tissueSiteDetailId,"&",
                 "datasetId=",datasetId
  )
  url1 <- utils::URLencode(url1)
  url1Get <- curl::curl_fetch_memory(url1)
  if(url1Get$status_code!=200){
    stop("Http status code: ", url1Get$status_code)
  }
  url1GetText <- rawToChar(url1Get$content)
  url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
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




