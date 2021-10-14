#' @title Query gene information
#' @description
#'  users can use this function to search gene of interest
#' @param genes Following gene types are supported:
#' \itemize{
#'   \item \strong{Gene symbol}.
#'
#'   A character string or a character vector (case ignored). like: "tp53","naDK","SDF4".
#'   \item \strong{Gencode/ensemble id} (versioned or unversioned).
#'
#'    A character string or a character vector (case ignored). like: "ENSG00000210195.2","ENSG00000078808".
#'   \item \strong{Entrez gene ID}.
#'
#'   A integer string or a integer vectors. like: 51150,5590,4509.
#' }
#'
#' @param geneType A character string. Types of queried genes. Options: "symbol" (default), "gencodeId", "entrezId";
#' @param gencodeVersion A character string. Two version are supported: "v26" (default) and "v19"
#' @param genomeBuild A character string. "GRCh38/hg38"(default) for gencodeVersion "v26", "GRCh37/hg19" for gencodeVersion "v19"
#' @import data.table
#' @import stringr
#' @return A data.table of queried gene information. including following columns:
#' \itemize{
#'  \item \strong{genes.} Input genes
#'  \item \strong{geneSymbol.} Gene symbol.
#'  \item \strong{gencodeId.} Gencode/ensemble id (versioned).
#'  \item \strong{entrezGeneId.} Entrez gene ID.
#'  \item \strong{geneType.} Gene type.
#'  \item \strong{chromosome.} Note: "chr" is added in gencode v26,
#'  \item \strong{start.}
#'  \item \strong{end.}
#'  \item \strong{strand.}
#'  \item \strong{tss.} Transcriptional start site.
#'  \item \strong{gencodeVersion.} Gencode Version.
#'  \item \strong{genomeBuild.} Genome version.
#'  \item \strong{description.}
#'  }
#'
#' @export
#'
#' @examples
#' \donttest{
#'  # test hg38:
#'  geneInfo <- GTExquery_gene("TP53", "symbol", "v26", "GRCh38/hg38")
#'  geneInfo <- GTExquery_gene(c("tp53","naDK","SDF4"), "symbol", "v26", "GRCh38/hg38")
#'  geneInfo <- GTExquery_gene(c("ENSG00000210195.2","ENSG00000078808"), geneType="gencodeId", "v26", "GRCh38/hg38")
#'  geneInfo <- GTExquery_gene(c(51150,5590,4509), "entrezId", "v26", "GRCh38/hg38")
#'  # test hg19:
#'  geneInfo <- GTExquery_gene(c(51150,5590,4509), "entrezId", "v19","GRCh37/hg19")
#'  }
GTExquery_gene <- function(genes="", geneType="symbol", gencodeVersion="v26", genomeBuild="GRCh38/hg38"){
  # check null/na
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }
  # geneType
  if( is.null(geneType) ||  is.na(geneType) || geneType=="" || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"symbol\", \"gencodeId\", \"entrezId\".")
  }else if( !(geneType %in% c("symbol", "gencodeId", "entrezId")) ){
    stop("Parameter \"geneType\" should be choosen from \"symbol\", \"gencodeId\", \"entrezId\".")
  }
  # gencodeVersion
  if( is.null(gencodeVersion) ||  is.na(gencodeVersion) || gencodeVersion=="" || length(gencodeVersion)!=1){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }else if( !(gencodeVersion %in% c("v26", "v19")) ){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }
  # genomeBuild
  # check gencodeVersion match with genomeBuild
  if( is.null(genomeBuild) ||  is.na(genomeBuild) || genomeBuild=="" || length(genomeBuild)!=1){
    stop("Parameter \"genomeBuild\" should be choosen from \"GRCh38/hg38\", \"GRCh37/hg19\".")
  }else if( !(genomeBuild %in% c("GRCh38/hg38", "GRCh37/hg19")) ){
    stop("Parameter \"genomeBuild\" should be choosen from \"GRCh38/hg38\", \"GRCh37/hg19\".")
  } else if(gencodeVersion=="v26" & genomeBuild=="GRCh38/hg38"){
    # data(gencodeGeneInfoV26)
    gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV26)
  }else if(gencodeVersion=="v19" & genomeBuild=="GRCh37/hg19"){
    # data(gencodeGeneInfoV19)
    gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV19)
  }else{
    stop("gencodeVersion must be matched with genomeBuild.\n eg. v26(GRCh38/hg38), v19(GRCh37/hg19)")
  }

  # merge:
  if( geneType=="symbol" ){
    genesDT <- data.table::data.table(genes=as.character(genes))
    genesDT$geneSymbolUpper <- toupper(genesDT$genes)
    gencodeGeneInfo$geneSymbolUpper<- toupper(gencodeGeneInfo$geneSymbol)
    genesDTout <- merge(genesDT, gencodeGeneInfo, by ="geneSymbolUpper", sort=FALSE)
    genesDTout <- genesDTout[,-c("geneSymbolUpper")]
    if( nrow(genesDTout) >0){
      if( length(genes)>1 ){
        message("Queried ",length(genes), " genes, finally ",nrow(genesDTout)," records matched!")
      }else if( length(genes)==1 ){
        message("Queried ",length(genes), " gene, finally ",nrow(genesDTout)," record matched!")
      }
      return(genesDTout)
    }else{
      if( length(genes)>1 ){
        message("Queried ",length(genes), " genes, finally 0 record matched!\nplease check your input genes!")
      }else if( length(genes)==1 ){
        message("Queried ",length(genes), " gene, finally 0 record matched!\nplease check your input genes!")
      }
      return(data.table::data.table(genes=genes))
    }
  }
  # gencodeId
  if( geneType=="gencodeId" ){
    genesDT <- data.table::data.table(genes=as.character(genes))
    genesDT$gencodeIdUpper <- toupper(genesDT$genes)
    gencodeGeneInfo$gencodeIdUpper<- toupper(gencodeGeneInfo$gencodeId)
    gencodeGeneInfo$gencodIdUv <- unlist(lapply(gencodeGeneInfo$gencodeId, function(x){ stringr::str_split_fixed(x, stringr::fixed("."),2)[1] }))
    gencodeGeneInfo$gencodIdUvUpper<- toupper(gencodeGeneInfo$gencodIdUv)
    genesDTout_1 <- merge(genesDT, gencodeGeneInfo, by ="gencodeIdUpper", sort=FALSE)
    genesDTout_2 <- merge( gencodeGeneInfo, genesDT, by.y ="gencodeIdUpper", by.x= "gencodIdUvUpper", sort=FALSE)
    genesDTout <- rbind(genesDTout_1, genesDTout_2)[,c("genes",names(gencodeGeneInfo[,-c("gencodeIdUpper", "gencodIdUv", "gencodIdUvUpper")])),with=FALSE]
    if( nrow(genesDTout) >0){
      message("Queried ",length(genes), " genes, finally ",nrow(genesDTout)," records matched!")
      return(genesDTout)
    }else{
      message("Queried ",length(genes), " genes, finally 0 records matched!\nplease check your input genes!")
      return(data.table::data.table(genes=genes))
    }
    # entrezId
  }else if( geneType=="entrezId" ){
    if( !all(unlist(lapply(genes, is.numeric))) ){
      stop("Integer is requeried for entrezId!")
    }
    genesDT <- data.table::data.table(genes=genes)
    gencodeGeneInfo$genes <- gencodeGeneInfo$entrezGeneId
    genesDTout <- merge(genesDT, gencodeGeneInfo, by ="genes", sort=FALSE)
    if( nrow(genesDTout) >0){
      message("Queried ",length(genes), " genes, finally ",nrow(genesDTout)," records matched!")
      return(genesDTout)
    }else{
      message("Queried ",length(genes), " genes, finally 0 records matched!\nplease check your input genes!")
      return(data.table::data.table(genes=genes))
    }
  }else{
    return(data.table::data.table(genes=genes))
  }
}

#' @title  Fetch information of samples used in analyses from all datasets.
#'
#' @param tissueSiteDetail
#' For GTEx v8, must be following terms:
#'  \itemize{
#'    \item All
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
#'   \item All
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
#' @param dataType A character string. Options: "RNASEQ" (default), "WGS", "WES", "OMNI".
#' @param datasetId A character string. Options: "gtex_v8" (default), "gtex_v7".
#' @param pageSize A integer value (1-2000). number of records fetched per request (default: 200).
#' @import data.table
#' @import utils
#' @import curl
#' @import jsonlite
#' @return returen sample information
#' @export
#' @examples
#' \donttest{
#'   GTExquery_sample( tissueSiteDetail="Liver", dataType="RNASEQ", datasetId="gtex_v8", pageSize=200 )
#'   GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ", datasetId="gtex_v8",pageSize=2000 )
#'   GTExquery_sample( "Adipose - Visceral (Omentum)", "RNASEQ", "gtex_v8", 200 )
#'   }
GTExquery_sample <- function( tissueSiteDetail="Liver", dataType="RNASEQ", datasetId="gtex_v8", pageSize=200 ){
  page_tmp <- 0
  pageSize_tmp <- as.integer(pageSize)

  ########## parameter check: tissueSiteDetail
  if(is.null(datasetId) ||  is.na(datasetId)){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(  !(datasetId %in% c("gtex_v8", "gtex_v7"))  ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }

  ########## parameter check: tissueSiteDetail
  # import tissueSiteDetailGTEx according to datasetId
  if(datasetId == "gtex_v8"){
    # data(tissueSiteDetailGTExv8)
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if(datasetId == "gtex_v7"){
    # data(tissueSiteDetailGTExv7)
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  is.na(tissueSiteDetail)){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if( !(tissueSiteDetail %in% c("All", tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("0. All", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }
  # convert tissueSiteDetail to tissueSiteDetailId:
  tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId


  ########## parameter check: pageSize_tmp
  if(is.null(pageSize_tmp) ||  is.na(pageSize_tmp)){
    stop("Parameter \"pageSize\" should be a integer!")
    }else if(length(pageSize_tmp)!=1){
      stop("Parameter \"pageSize\" should be a integer!")
    }else if( pageSize_tmp > 2000 | pageSize_tmp < 1){
      stop("pageSize must between 1 and 2000!")
    }

  ######## parameter check: dataType
  if(is.null(dataType) ||  is.na(dataType)){
    stop("Parameter \"dataType\" should be a character!")
  }else if(length(dataType)!=1){
    stop("Parameter \"dataType\" should be a character!")
  }else if( !dataType %in% c("RNASEQ", "WGS", "WES", "OMNI")){ # "EXCLUDE", "USE ME"
    stop("Parameter \"dataType\" should be chosen from: \"RNASEQ\", \"WGS\", \"WES\", \"OMNI\"." )
  }


  ########## construct api url:
  if( tissueSiteDetail =="All"){
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                   "datasetId=", datasetId,"&",
                   "tissueSiteDetail=", "All","&",
                   "dataType=",dataType,"&",
                   "page=",page_tmp,"&",
                   "pageSize=", pageSize_tmp, "&",
                   "sortBy=sampleId&sortDirection=asc"
    )
  }else{
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                   "datasetId=", datasetId,"&",
                   "tissueSiteDetailId=", tissueSiteDetailId,"&",
                   "dataType=",dataType,"&",
                   "page=",page_tmp,"&",
                   "pageSize=", pageSize_tmp, "&",
                   "sortBy=sampleId&sortDirection=asc"
    )
  }
  url1 <- utils::URLencode(url1)
  # url1 <- "https://gtexportal.org/rest/v1/dataset/sample?datasetId=gtex_v8&tissueSiteDetailId=Liver&dataType=RNASEQ&format=json&page=0&pageSize=2000&sortBy=sampleId&sortDirection=asc"
  # url1 <- "https://gtexportal.org/rest/v1/dataset/sample?datasetId=gtex_v8&tissueSiteDetail=All&dataType=RNASEQ&format=json&page=0&pageSize=200&sortBy=sampleId&sortDirection=asc"
  outInfo <- data.table::data.table()
  pingOut <- apiAdmin_ping()
  if( !is.null(pingOut) && pingOut==200 ){
    message("GTEx API successfully accessed!")
    url1Get <- curl::curl_fetch_memory(url1)
    url1GetText <- rawToChar(url1Get$content)
    url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
    tmp <- data.table::as.data.table(url1GetText2Json$sample)
    outInfo <- rbind(outInfo, tmp)
    message("Total records: ",url1GetText2Json$recordsFiltered,"; total pages: ",url1GetText2Json$numPages,"; downloaded pages:",url1GetText2Json$page," records: ", nrow(tmp))
    while(url1GetText2Json$page < (url1GetText2Json$numPages-1)){
      page_tmp <- page_tmp+1
      if( tissueSiteDetail =="All"){
        url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                       "datasetId=", datasetId,"&",
                       "tissueSiteDetail=", "All","&",
                       "dataType=",dataType,"&",
                       "page=",page_tmp,"&",
                       "pageSize=", pageSize_tmp, "&",
                       "sortBy=sampleId&sortDirection=asc"
        )
      }else{
        url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                       "datasetId=", datasetId,"&",
                       "tissueSiteDetailId=", tissueSiteDetailId,"&",
                       "dataType=",dataType,"&",
                       "page=",page_tmp,"&",
                       "pageSize=", pageSize_tmp, "&",
                       "sortBy=sampleId&sortDirection=asc"
        )
      }
      url1 <- utils::URLencode(url1)
      url1Get <- curl::curl_fetch_memory(url1)
      url1GetText <- rawToChar(url1Get$content)
      url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
      tmp <- data.table::as.data.table(url1GetText2Json$sample)
      outInfo <- rbind(outInfo, tmp)
      message("Total records: ",url1GetText2Json$recordsFiltered,"; total pages: ",url1GetText2Json$numPages,"; downloaded pages:",url1GetText2Json$page," records: ", nrow(tmp))
      return(outInfo[,.(sampleId, sex, ageBracket,datasetId, tissueSiteDetail, tissueSiteDetailId, pathologyNotes, hardyScale,  dataType )])
    }
  }else{
    message("")
  }
  # note: pathologyNotes info is ignored.
}


#' @title Fetch normalized gene expression.
#' @description
#'  This function fetch queried genes' expression profiles in a tissue at the sample level.
#' @param genes Following gene types are supported:
#' \itemize{
#'   \item \strong{Gene symbol}.
#'
#'   A character string or a character vector (case ignored). like: "tp53","naDK","SDF4".
#'   \item \strong{Gencode/ensemble id} (versioned or unversioned).
#'
#'    A character string or a character vector (case ignored). like: "ENSG00000210195.2","ENSG00000078808".
#'   \item \strong{Entrez gene ID}.
#'
#'   A integer string or a integer vectors. like: 51150,5590,4509.
#' }
#'
#' @param geneType A character string. Types of queried genes. Options: "symbol" (default), "gencodeId", "entrezId";
#' @param tissueSiteDetail
#'  For GTEx v8, must be following terms:
#'  \itemize{
#'    \item All
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
#'   \item All
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
#' @import data.table
#' @import curl
#' @import stringr
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom   IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @return return gene expression profiles.
#' @export
#' @examples
#' \donttest{
#'   expProfiles <- GTExquery_exp(c("ENSG00000210195.2","ENSG00000078808"), "gencodeId", "Liver", "gtex_v8")
#'   expProfiles <- GTExquery_exp(c("tp53","naDK","SDF4"), "symbol", "Artery - Coronary", "gtex_v8", toSummarizedExperiment=TRUE)
#'   expProfiles <- GTExquery_exp(c("tp53","naDK","SDF4"), "symbol", "Artery - Coronary", "gtex_v7")
#'   }
GTExquery_exp <- function(genes="", geneType="symbol", tissueSiteDetail="Liver", datasetId="gtex_v8", toSummarizedExperiment=FALSE ){
  # genes=c("ENSG00000210195.2","ENSG00000078808")
  # geneType="gencodeId"
  # tissueSiteDetail="Liver"
  # datasetId="gtex_v8"
  # toSummarizedExperiment=TRUE

  ########## parameter check: tissueSiteDetail
  if(is.null(datasetId) ||  is.na(datasetId)){
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
    gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV26)
    gencodeVersion <- "v26"
    genomeBuild <- "GRCh38/hg38"
  }else if(datasetId == "gtex_v7"){
    # data(tissueSiteDetailGTExv7)
    # data(gencodeGeneInfoV19)
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
    gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV19)
    gencodeVersion <- "v19"
    genomeBuild <- "GRCh37/hg19"
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  is.na(tissueSiteDetail)){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if( !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }
  # convert tissueSiteDetail to tissueSiteDetailId:
  tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId

  # convert genes. parameter check is unnecessary for this, because GTExquery_gene check it internally.
  message("Querying gene info:")
  geneInfo <- GTExquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, genomeBuild)
  # geneInfo <-  GTExquery_gene(c("tp53","naDK","SDF4"), "symbol", "v26", "GRCh38/hg38")

  # get sample info:
  message("Fetching sample information from API server:")
  sampleInfo <- GTExquery_sample(tissueSiteDetail=tissueSiteDetail, dataType="RNASEQ",datasetId=datasetId, pageSize=200 )
  if(!exists("sampleInfo")){
    stop("can not connect to API server, please try again later.")
  }


  ################# construct url:
  if( tissueSiteDetail =="All"){
    url1 <- paste0("https://gtexportal.org/rest/v1/expression/geneExpression?",
                   "datasetId=", datasetId,"&",
                   "gencodeId=", paste(geneInfo$gencodeId, collapse = ","), "&",
                   "tissueSiteDetail=", "All","&",
                   "&format=json"
    )
  }else{
    url1 <- paste0("https://gtexportal.org/rest/v1/expression/geneExpression?",
                   "datasetId=", datasetId,"&",
                   "gencodeId=", paste(geneInfo$gencodeId, collapse = ","), "&",
                   "tissueSiteDetailId=", tissueSiteDetailId,"&",
                   "&format=json"
    )
  }
  url1 <- utils::URLencode(url1)
  # url1 <- "https://gtexportal.org/rest/v1/expression/geneExpression?datasetId=gtex_v8&gencodeId=ENSG00000240361.1%2CENSG00000222623.1%2CENSG00000008128.22%2CENSG00000008130.15&tissueSiteDetailId=Bladder&format=json"
  # url1 <- "https://gtexportal.org/rest/v1/expression/geneExpression?datasetId=gtex_v8&gencodeId=ENSG00000240361.1%2CENSG00000222623.1%2CENSG00000008128.22%2CENSG00000008130.15&tissueSiteDetail=All&format=json"
  outInfo <- data.table::data.table()
  pingOut <- apiAdmin_ping()
  if( !is.null(pingOut) && pingOut==200 ){
    message("GTEx API successfully accessed!")
    url1Get <- curl::curl_fetch_memory(url1)
    url1GetText <- rawToChar(url1Get$content)
    url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
    outInfo <- data.table::as.data.table(url1GetText2Json$geneExpression)
    outInfo <- merge(geneInfo, outInfo, by= c("geneSymbol", "gencodeId"))
    outInfo <- cbind(outInfo[,c("genes")], outInfo[,-c("genes")])
    expDT <- data.table::as.data.table(t(data.table::as.data.table(outInfo$data)))
    data.table::setDF(expDT)
    rownames(expDT) <- geneInfo$geneSymbol
    colnames(expDT) <- sampleInfo$sampleId
    message( length(outInfo$data[[1]])," Expression profiles of ",length(unique(outInfo$gencodeId)), " genes in ",unique(outInfo$tissueSiteDetailId)," were obtained."," Unit: ",paste0(unique(outInfo$unit), collapse = ","),"." )
    #
    if( !toSummarizedExperiment ){
      outInfoDF <- cbind(outInfo[,-c("data")], expDT)
      data.table::setDF(outInfoDF)
      return(outInfoDF)
    }else{
      # construc summarizedExperiment object::
      # gene info:
      expRowRanges <- GenomicRanges::GRanges(geneInfo$chromosome,
                                          IRanges::IRanges(geneInfo$start, geneInfo$end),
                                          strand= geneInfo$strand,
                                          geneInfo[,.(genes, geneSymbol, gencodeId, entrezGeneId, geneType, tss, description)]
                                          # ,seqinfo = GenomeInfoDb::Seqinfo(genome= stringr::str_split(unique(geneInfo$genomeBuild)[1],stringr::fixed("/"))[[1]][2])
                                          )
      # sample info:
      expColData <- data.table::copy(sampleInfo)
      # meta data:
      expMetadata <- paste0(length(outInfo$data[[1]])," Expression profiles of ",length(unique(outInfo$gencodeId)), " genes in ",unique(outInfo$tissueSiteDetailId), ". Unit: ",paste0(unique(outInfo$unit), collapse = ","),".")
      # outInfoSE:
      outInfoSE <- SummarizedExperiment::SummarizedExperiment(assays=list(exp=expDT),
                                        rowRanges=expRowRanges,
                                        colData=expColData,
                                        metadata=expMetadata)
    }

  }else{
    message("")
  }
}


#
#' @title Heartbeat to check server connectivity.
#' @description
#'  test API server
#' @import curl
#' @return A boolean value.
#' @export
#' @examples
#' \donttest{
#'   apiAdmin_ping()
#'   apiStatus <- ifelse( apiAdmin_ping() ==200, "accessed", "failed")
#'   print(apiStatus)
#'  }
apiAdmin_ping <- function(){
  url1Get <- "https://gtexportal.org/rest/v1/admin/ping"
  tryCatch(
    {
      # httr::status_code(httr::GET(url1Get))
      curl::curl_fetch_memory(url1Get)$status_code
    },
    # e = simpleError("test error"),
    error=function(cond){
      message("Note: API server is busy or your network has latency, please try again later.")
      return(NULL)
    },
    warning = function(cond){
      message(cond)
    }
  )
}














