
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
#'                                  "gencodeId", "Liver", "gtex_v8")
#'   # Download gene expression profiles with multiple genes:
#'   expProfiles <- xQTLdownload_exp(c("tp53","naDK","SDF4"),
#'                                  "geneSymbol", "Artery - Coronary", "gtex_v8",
#'                                  pathologyNotesCategories=TRUE,
#'                                  toSummarizedExperiment=TRUE)
#'   expProfiles <- xQTLdownload_exp(c("tp53","naDK","SDF4"),
#'                                  "geneSymbol", "Artery - Coronary", "gtex_v7")
#'
#'   # Get proteing-coding genes' expression:
#'   proT <- xQTLquery_gene (genes="protein coding", geneType="geneCategory", "v26" )
#'   proTexp <- xQTLdownload_exp(proT$geneSymbol[1:100], geneType = "geneSymbol","Lung","gtex_v8")
#'   }
xQTLdownload_exp <- function(genes="", geneType="geneSymbol", tissueSiteDetail="Liver", datasetId="gtex_v8", toSummarizedExperiment=TRUE, recordPerChunk=150, pathologyNotesCategories=FALSE  ){
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

  ############ convert genes. parameter check is unnecessary for this, because xQTLquery_gene check it internally.
  message("== Querying gene info from API server:")
  geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, recordPerChunk=recordPerChunk)
  if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
    stop("gene information is null.")
  }
  message("== Done.")
  # geneInfo <-  xQTLquery_gene(c("tp53","naDK","SDF4"), "geneSymbol", "v26", "GRCh38/hg38")

  ############ get sample info:
  message("== Fetching sample information from API server:")
  sampleInfo <- xQTLquery_sample(tissueSiteDetail=tissueSiteDetail, dataType="RNASEQ", datasetId=datasetId, recordPerChunk=recordPerChunk,pathologyNotesCategories=pathologyNotesCategories )
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
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
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
    url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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
#'  xQTLdownload_eqtlSig(variantName="rs201327123", variantType="snpId")
#'  xQTLdownload_eqtlSig(variantName="chr1_14677_G_A_b38", variantType="variantId")
#'  xQTLdownload_eqtlSig(variantName="11_66328719_T_C_b37", variantType="variantId",
#'                    datasetId="gtex_v7")
#'  xQTLdownload_eqtlSig(variantName="11_66328719_T_C_b37", variantType="variantId",
#'                    datasetId="gtex_v7", tissueSiteDetail="Skin - Sun Exposed (Lower leg)")
#'
#'  # Download eQTL association according to all tissues with genome location:
#'  varInfo <-  xQTLquery_varPos(chrom="chr1", pos=c(1102708),"gtex_v8")
#'  xQTLdownload_eqtlSig(variantName=varInfo$snpId, variantType="snpId")
#'
#'  # Download eQTL info for a gene:
#'  xQTLdownload_eqtlSig(gene="ATAD3B")
#'  xQTLdownload_eqtlSig(gene="TP53",datasetId="gtex_v7")
#'  xQTLdownload_eqtlSig(gene="ENSG00000141510.16", geneType="gencodeId", datasetId="gtex_v8")
#'  xQTLdownload_eqtlSig(gene="ENSG00000141510.11", geneType="gencodeId",
#'                    datasetId="gtex_v7",tissueSiteDetail="Thyroid" )
#'  xQTLdownload_eqtl(gene="ENSG00000141510.11", geneType="gencodeId",
#'                    datasetId="gtex_v7",tissueSiteDetail="Thyroid" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  xQTLdownload_eqtlSig(variantName="rs1641513",gene="TP53", datasetId="gtex_v8")
#'  xQTLdownload_eqtlSig(variantName="rs1641513",gene="TP53", datasetId="gtex_v7")
#'  xQTLdownload_eqtlSig(variantName="chr1_1667948_A_G_b38", variantType="variantId",
#'                    gene="SLC35E2B", geneType="geneSymbol", tissueSiteDetail="Kidney - Cortex")
#' }
xQTLdownload_eqtlSig <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8"){
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

  ##################### fetch varInfo
  if(variantName!=""){
    message("== Querying variant info from API server:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      message("The variant [",variantName, "] is not incuded in [",datasetId,"].")
      return(NULL)
    }else{
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if(gene !=""){
    message("== Querying gene info from API server:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input.")
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
  # check network:
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  tmp <- data.table::as.data.table(url1GetText2Json$singleTissueEqtl)
  if( !exists("tmp")||nrow(tmp)==0){
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
#' @param recordPerChunk A integer value (1-500). number of records fetched per request (default: 100).
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
#'  a <- xQTLdownload_eqtl(gene="TP53")
#'  xQTLdownload_eqtl(gene="ATAD3B", datasetId="gtex_v7")
#'
#'  # Unversioned gencode ID in GTEx V7:
#'  eqtl_v7 <- xQTLdownload_eqtl(gene="ENSG00000141510", geneType="gencodeId",
#'                                  datasetId="gtex_v7")
#'  # Unversioned gencode ID in GTEx V8:
#'  eqtl_v8 <- xQTLdownload_eqtl(gene="ENSG00000141510", geneType="gencodeId",
#'                                  datasetId="gtex_v8")
#'
#'  # In a specific tissue:
#'  xQTLdownload_eqtl(gene="ENSG00000141510.11", geneType="gencodeId",
#'                       datasetId="gtex_v7", tissueSiteDetail="Thyroid" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  xQTLdownload_eqtl(variantName="rs1641513",gene="TP53", datasetId="gtex_v8")
#'  xQTLdownload_eqtl(variantName="rs11657498",gene="TP53",
#'                       datasetId="gtex_v7")
#'  xQTLdownload_eqtlSig(variantName="chr1_1667948_A_G_b38", variantType="variantId",
#'                       gene="SLC35E2B", geneType="geneSymbol",
#'                       tissueSiteDetail="Kidney - Cortex")
#'  xQTLdownload_eqtl(variantName="17_7492388_G_A_b37", variantType="variantId",
#'                       gene="TP53", geneType="geneSymbol",
#'                       tissueSiteDetail="Uterus",  datasetId="gtex_v7")
#' }
xQTLdownload_eqtl <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8", recordPerChunk=100){
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
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  # message("GTEx API successfully accessed!")

  ##################### fetch geneInfo:
  if(gene !=""){
    message("== Querying gene info from API server:")
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
    message("== Querying variant info from API server:")
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
  url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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

#' @title Fetch all eQTL associations from EBI eQTL category.
#'
#' @param gene gene
#' @param geneType geneType
#' @param tissueSiteDetail tissueSiteDetail
#' @param recordPerChunk A integer value (1-200). number of records fetched per request (default: 200).
#' @param study "GTEx_V8"
#' @param withB37VariantId Whether to return the genome location(GTEx v7) of variants. Default: TRUE.
#' @import data.table
#' @import stringr
#' @return a data.table
#' @export
#'
#' @examples
#' \donttest{
#'   geneAsso <- xQTLdownload_assoAll("ATP11B", tissueSiteDetail="Muscle - Skeletal", withB37VariantId=FALSE)
#' }
xQTLdownload_assoAll <- function(gene="", geneType="geneSymbol", tissueSiteDetail="", recordPerChunk=250, study="GTEx_V8", withB37VariantId=TRUE){
  .<-NULL
  variantId <- variant <- b37VariantId <- snpId <- NULL
  # gene="CYP2W1"
  # geneType="geneSymbol"
  # tissueSiteDetail="Lung"

  study="GTEx_V8"

  ########## check genes
  if( is.null(gene) ||  any(is.na(gene)) || any(gene=="") ||length(gene)==0 ){
    stop("Parameter \"gene\" can not be NULL or NA!")
  }else if( length(gene)!=1 ){
    stop("Gene list is not supported, please input a gene symbol or gencode ID.")
  }

  ########## geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("geneSymbol", "gencodeId", "geneCategory")) ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
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
    message("== Querying gene info from API server:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion)
    geneInfoV19 <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v19")
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or leave \"gene\" as null")
    }else if( nrow(geneInfo)>1 || nrow(geneInfoV19)>1 ){
      stop("Totally, ",nrow(geneInfo), " gencode ID of queried gene [", gene,"] were detected, please enter the gencode ID (versioned or unversioned) for querying!")
    }else{
      message("== Done.")
    }
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
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'   xQTLdownload_aQTLAllAsso( gene="TP53", tissueSiteDetail = "Adipose - Visceral (Omentum)" )
#'   xQTLdownload_aQTLAllAsso( refSeq = "NM_001291581.2", tissueSiteDetail = "Adipose - Visceral (Omentum)" )
#'   xQTLdownload_aQTLAllAsso( rsid = "rs3026133", tissueSiteDetail = "Adipose - Visceral (Omentum)", pThreshold=1e-1 )
#'   xQTLdownload_aQTLAllAsso( variantId = "chr17_5289661_A_G_b38", tissueSiteDetail = "Adipose - Visceral (Omentum)" )
#' }
xQTLdownload_aQTLAllAsso<- function(gene="", geneType="geneSymbol", refSeq="", rsid="",  variantId="", pThreshold=1, tissueSiteDetail=""){
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
#'  xQTLdownload_sqtlSig(variantName="rs201327123", variantType="snpId")
#'  xQTLdownload_sqtlSig(variantName="chr1_14677_G_A_b38", variantType="variantId")
#'  xQTLdownload_sqtlSig(variantName="chr1_14677_G_A_b38", variantType="variantId",
#'                       tissueSiteDetail="Whole Blood")
#'
#'  # Download sQTL association according to all tissues with genome location:
#'  varInfo <-  xQTLquery_varPos(chrom="chr1", pos=c(1102708),"gtex_v8")
#'  xQTLdownload_sqtlSig(variantName=varInfo$snpId, variantType="snpId")
#'
#'  # Download sQTL info for a gene:
#'  xQTLdownload_sqtlSig(gene="ATAD3B")
#'  xQTLdownload_sqtlSig(gene="ENSG00000141510.16", geneType="gencodeId")
#'  xQTLdownload_sqtlSig(gene="ENSG00000141510.16", geneType="gencodeId",
#'                       tissueSiteDetail="Lung" )
#'
#'  # Download eQTL info for a variant-gene pair:
#'  xQTLdownload_sqtlSig(variantName="rs546057177", gene="TP53")
#'  xQTLdownload_sqtlSig(variantName="chr17_7465085_A_G_b38", variantType="variantId",
#'                       gene="TP53", geneType="geneSymbol",
#'                       tissueSiteDetail="Lung")
#' }
xQTLdownload_sqtlSig <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="" ){
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
    message("== Querying variant info from API server:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType, datasetId=datasetId)
    if(nrow(varInfo)==0 || is.null(varInfo)|| !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }

  ##################### fetch geneInfo:
  if(gene !=""){
    message("== Querying gene info from API server:")
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
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  url1 <- utils::URLencode(url1)
  url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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
#'                       variantType="variantId", geneType="gencodeId",
#'                      tissueSiteDetail="Stomach")
#' }
xQTLdownload_eqtlExp <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8"){
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
    message("== Querying variant info from API server:")
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
    message("== Querying gene info from API server:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
    if(nrow(geneInfo)==0|| is.null(geneInfo) || !exists("geneInfo")){
      stop("Invalid gene name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"gene\" can not be null!")
  }

  message("== Querying expression from API server:")
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  # construct url
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/rest/v1/association/dyneqtl?",
                 "gencodeId=",geneInfo$gencodeId,"&",
                 "variantId=",varInfo$variantId,"&",
                 "tissueSiteDetailId=",tissueSiteDetailId,"&",
                 "datasetId=",datasetId
  )
  url1 <- utils::URLencode(url1)
  url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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


#' @title Fetch linkage disequilibrium data
#' @description
#'   Fetch linkage disequilibrium data for the cis-eQTLs found associated with this gene in a specified dataset.
#' @param gene A gene symbol, gencode id (versioned), or a charater string of gene type.
#' @param geneType  A character string. "geneSymbol"(default), or "gencodeId".
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
#'  aa <- xQTLdownload_ld("TP53" )
#'  xQTLdownload_ld(gene = "DDX11", datasetId="gtex_v7")
#'  xQTLdownload_ld(gene="ENSG00000008128.22", geneType="gencodeId")
#'
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
  message("== Querying gene info from API server:")
  geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType, gencodeVersion=gencodeVersion)
  if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
    stop("The gene [",gene,"] you entered could not be found!")
  }
  message("== Done.")

  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  url1 <- paste0("https://gtexportal.org/rest/v1/dataset/ld?format=json&",
                 "gencodeId=",geneInfo$gencodeId,
                 "&datasetId=", datasetId)
  url1 <- utils::URLencode(url1)
  message("== Downloading...")
  url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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


#' @title Fetch eGenes (eQTL Genes) from the specified dataset.
#' @description
#'  eGenes are genes that have at least one significant cis-eQTL acting upon them. Results may be filtered by tissue.
#' @param gene  A charater string of gene symbol, gencode id (versioned).
#' @param geneType A character string. "geneSymbol"(default) or"gencodeId".
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
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
#' @param recordPerChunk A integer value (1-200). number of records fetched per request (default: 200).
#' @import data.table
#' @import stringr
#' @import utils
#' @return a data.table
#' @export
#'
#' @examples
#' \donttest{
#'  eGeneInfo <- xQTLdownload_egene("TP53")
#'  eGeneInfo <- xQTLdownload_egene(tissueSiteDetail="Lung", recordPerChunk=2000)
#'  eGeneInfo <- xQTLdownload_egene("ENSG00000141510.16", geneType="gencodeId")
#'  eGeneInfo <- xQTLdownload_egene("DDX11", datasetId="gtex_v7", tissueSiteDetail="Artery - Tibial" )
#' }
xQTLdownload_egene <- function(gene = "", geneType="geneSymbol", datasetId = "gtex_v8", tissueSiteDetail="", recordPerChunk=200){
  .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss <- log2AllelicFoldChange <- empiricalPValue <- pValue <- pValueThreshold <- qValue <-NULL
  # gene="DDX11"
  # geneType="geneSymbol"
  # datasetId = "gtex_v8"
  # tissueSiteDetail="Lung"
  # recordPerChunk=100

  page_tmp <- 0
  pageSize_tmp <- recordPerChunk

  # check genes
  if( is.null(gene) ||  any(is.na(gene)) || length(gene)==0 ){
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
    message("== Querying gene info from API server:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType, gencodeVersion=gencodeVersion)
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",gene,"] you entered could not be found!")
    }
    message("== Done.")
  }

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
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  # message("GTEx API successfully accessed!")
  suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
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
    suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
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


#' @title fetch gene median expression
#'
#' @param genes A character vector of gene symbol/gencode ID.
#' @param geneType "geneSymbol" or "gencodeID".
#' @param datasetId "gtex_v7" or "gtex_v8"
#' @param tissueSiteDetail tissue name.
#' @param recordPerChunk 1-2000. Default: 150.
#' @import data.table
#' @import utils
#' @import stringr
#' @return A data.frame
#' @export
#'
#' @examples
#' \donttest{
#'  xQTLdownload_geneMedExp(genes="TP53", "geneSymbol")
#' }
xQTLdownload_geneMedExp <- function(genes="", geneType="geneSymbol", datasetId="gtex_v8", tissueSiteDetail="", recordPerChunk=150 ){
  .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss<-strand <- NULL
  # check genes
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }

  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
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
  if(genes!="" || length(genes)>1){
    message("== Querying gene info from API server:")
    suppressMessages( geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, recordPerChunk = recordPerChunk))
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",stringr::str_sub(paste(genes, collapse = ","),1,20),"....] you entered could not be found!")
    }
    message("== Done.")
  }

  #
  outInfo <- data.table::data.table()
  genesCut <- data.table::data.table(gencodeId=geneInfo$gencodeId, ID=1:nrow(geneInfo), cutF = as.character(cut(1:nrow(geneInfo),breaks=seq(0,nrow(geneInfo)+recordPerChunk,recordPerChunk) )) )
  genesURL <- genesCut[,.(genesURL=paste0(gencodeId,collapse = "%2C")),by=c("cutF")]
  if( any(unlist(lapply(genesURL$genesURL, nchar)) >3900) ){
    stop("Too many queried genes, please lower the value of \"recordPerChunk\", or reduce your input genes.")
  }
  outInfo <- data.table::data.table()
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
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
    url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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
#'  retrieveLD('6', 'rs9349379', 'AFR')
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
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'   snpLD <- retrieveLD_ldlink("rs3", windowSize=50000)
#' }
retrieveLD_ldlink <- function(targetSnp="", population="EUR" , windowSize=500000, method="download", genomeVersion="grch38", max_count=3, token="9246d2db7917"){
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
