
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
#'  For GTEx v8, must be chosen from the following terms:
#'  \itemize{
#'    \item Adipose - Subcutaneous
#'    \item Adipose - Visceral (Omentum)
#'    \item Adrenal Gland
#'    \item Artery - Aorta
#'    \item Artery - Coronary
#'    \item Artery - Tibial
#'    \item Bladder
#'    \item Brain - Amygdala
#'    \item Brain - Anterior cingulate cortex (BA24)
#'    \item Brain - Caudate (basal ganglia)
#'    \item Brain - Cerebellar Hemisphere
#'    \item Brain - Cerebellum
#'    \item Brain - Cortex
#'    \item Brain - Frontal Cortex (BA9)
#'    \item Brain - Hippocampus
#'    \item Brain - Hypothalamus
#'    \item Brain - Nucleus accumbens (basal ganglia)
#'    \item Brain - Putamen (basal ganglia)
#'    \item Brain - Spinal cord (cervical c-1)
#'    \item Brain - Substantia nigra
#'    \item Breast - Mammary Tissue
#'    \item Cells - Cultured fibroblasts
#'    \item Cells - EBV-transformed lymphocytes
#'    \item Cervix - Ectocervix
#'    \item Cervix - Endocervix
#'    \item Colon - Sigmoid
#'    \item Colon - Transverse
#'    \item Esophagus - Gastroesophageal Junction
#'    \item Esophagus - Mucosa
#'    \item Esophagus - Muscularis
#'    \item Fallopian Tube
#'    \item Heart - Atrial Appendage
#'    \item Heart - Left Ventricle
#'    \item Kidney - Cortex
#'    \item Kidney - Medulla
#'    \item Liver
#'    \item Lung
#'    \item Minor Salivary Gland
#'    \item Muscle - Skeletal
#'    \item Nerve - Tibial
#'    \item Ovary
#'    \item Pancreas
#'    \item Pituitary
#'    \item Prostate
#'    \item Skin - Not Sun Exposed (Suprapubic)
#'    \item Skin - Sun Exposed (Lower leg)
#'    \item Small Intestine - Terminal Ileum
#'    \item Spleen
#'    \item Stomach
#'    \item Testis
#'    \item Thyroid
#'    \item Uterus
#'    \item Vagina
#'    \item Whole Blood
#'  }
#'  For GTEx v7, must be following terms:
#'  \itemize{
#'   \item Adipose - Subcutaneous
#'   \item Adipose - Visceral (Omentum)
#'   \item Adrenal Gland
#'   \item Artery - Aorta
#'   \item Artery - Coronary
#'   \item Artery - Tibial
#'   \item Bladder
#'   \item Brain - Amygdala
#'   \item Brain - Anterior cingulate cortex (BA24)
#'   \item Brain - Caudate (basal ganglia)
#'   \item Brain - Cerebellar Hemisphere
#'   \item Brain - Cerebellum
#'   \item Brain - Cortex
#'   \item Brain - Frontal Cortex (BA9)
#'   \item Brain - Hippocampus
#'   \item Brain - Hypothalamus
#'   \item Brain - Nucleus accumbens (basal ganglia)
#'   \item Brain - Putamen (basal ganglia)
#'   \item Brain - Spinal cord (cervical c-1)
#'   \item Brain - Substantia nigra
#'   \item Breast - Mammary Tissue
#'   \item Cells - EBV-transformed lymphocytes
#'   \item Cells - Transformed fibroblasts
#'   \item Cervix - Ectocervix
#'   \item Cervix - Endocervix
#'   \item Colon - Sigmoid
#'   \item Colon - Transverse
#'   \item Esophagus - Gastroesophageal Junction
#'   \item Esophagus - Mucosa
#'   \item Esophagus - Muscularis
#'   \item Fallopian Tube
#'   \item Heart - Atrial Appendage
#'   \item Heart - Left Ventricle
#'   \item Kidney - Cortex
#'   \item Liver
#'   \item Lung
#'   \item Minor Salivary Gland
#'   \item Muscle - Skeletal
#'   \item Nerve - Tibial
#'   \item Ovary
#'   \item Pancreas
#'   \item Pituitary
#'   \item Prostate
#'   \item Skin - Not Sun Exposed (Suprapubic)
#'   \item Skin - Sun Exposed (Lower leg)
#'   \item Small Intestine - Terminal Ileum
#'   \item Spleen
#'   \item Stomach
#'   \item Testis
#'   \item Thyroid
#'   \item Uterus
#'   \item Vagina
#'   \item Whole Blood
#'  }
#' @param datasetId A character string. Options: "gtex_v8" (default), "gtex_v7".
#' @param toSummarizedExperiment whether to return a data.frame or a summarizedExperiment object. Default: TRUE, return a toSummarizedExperiment object.
#' @param recordPerChunk A integer value (1-2000). number of records fetched per request (default: 150).
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom   IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @return return a SummarizedExperiment or data.frame object harbing gene expression profiles.
#' @export
#' @examples
#' \donttest{
#'   expProfiles <- GTExdownload_exp(c("ENSG00000210195.2"),
#'                                  "gencodeId", "Liver", "gtex_v8")
#'   expProfiles <- GTExdownload_exp(c("tp53","naDK","SDF4"),
#'                                  "geneSymbol", "Artery - Coronary", "gtex_v8",
#'                                  toSummarizedExperiment=TRUE)
#'   expProfiles <- GTExdownload_exp(c("tp53","naDK","SDF4"),
#'                                  "geneSymbol", "Artery - Coronary", "gtex_v7")
#'   # get miRNA expression:
#'   miRNA <- GTExquery_gene (genes="miRNA", geneType="geneCategory" )
#'   miRNAExp <- GTExdownload_exp(miRNA$gencodeId, geneType = "gencodeId")
#'   geneExp <- GTExdownload_exp(c("tp53",miRNA$geneSymbol[1:100],"naDK"), geneType = "geneSymbol")
#'   proT <- GTExquery_gene (genes="protein coding", geneType="geneCategory", "v26" )
#'   proTexp <- GTExdownload_exp(proT$geneSymbol[1:100], geneType = "geneSymbol","Lung","gtex_v7")
#'   }
GTExdownload_exp <- function(genes="", geneType="geneSymbol", tissueSiteDetail="Liver", datasetId="gtex_v8", toSummarizedExperiment=TRUE, recordPerChunk=150  ){
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
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
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
