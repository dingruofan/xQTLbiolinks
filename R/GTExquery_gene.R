#' @title query gene info
#' @description
#'  users can use this function to search gene of interest
#' @param genes gene names or gene ids (case ignored). Following gene types are supported: gene symbol (character); gencode/ensemble id(versioned or unversioned)(character); entrez gene ID (integer); Default: gene symbol (like: "tp53","naDK","SDF4").
#' @param geneType types of queried genes. Options: "symbol", "gencodeId", "entrezId";
#' @param gencodeVersion two version are supported: "v26" and "v19"
#' @param genomeBuild "GRCh38/hg38" for gencodeVersion "v26", "GRCh37/hg19" for gencodeVersion "v19"
#' @import data.table
#' @import stringr
#' @return
#' @export
#'
#' @examples
#'  # test hg38:
#'  GTExquery_gene("TP53", "symbol", "v26", "GRCh38/hg38")
#'  GTExquery_gene(c("tp53","naDK","SDF4"), "symbol", "v26", "GRCh38/hg38")
#'  GTExquery_gene(c("ENSG00000210195.2","ENSG00000078808","Ensg00000008130.15"), "symbol", "v26", "GRCh38/hg38")
#'  GTExquery_gene(c(51150,5590,4509), "entrezId", "v26", "GRCh38/hg38")
#'  # test hg19:
#'  GTExquery_gene(c(51150,5590,4509), "entrezId", "v19","GRCh37/hg19")
GTExquery_gene <- function(genes, geneType="symbol", gencodeVersion="v26", genomeBuild="GRCh38/hg38"){
  genes <- unlist(genes)
  if( any(is.na(genes)) | any(genes=="") | any(is.null(genes)) |length(genes)==0 ){
    stop("Gene name/ID can not be null!")
  }else if(gencodeVersion=="v26" & genomeBuild=="GRCh38/hg38"){
    # data(gencodeGeneInfoV26)
    gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV26)
  }else if(gencodeVersion=="v19" & genomeBuild=="GRCh37/hg19"){
    # data(gencodeGeneInfoV19)
    gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV19)
  }else{
    stop("gencodeVersion must be matched with genomeBuild.\neg. v26(GRCh38/hg38), v19(GRCh37/hg19)")
  }

  # merge:
  if( geneType=="symbol" ){
    genesDT <- data.table::data.table(genes=as.character(genes))
    genesDT$geneSymbolUpper <- toupper(genesDT$genes)
    gencodeGeneInfo$geneSymbolUpper<- toupper(gencodeGeneInfo$geneSymbol)
    genesDTout <- merge(genesDT, gencodeGeneInfo, by ="geneSymbolUpper", sort=FALSE)
    genesDTout <- genesDTout[,-c("geneSymbolUpper")]
    if( nrow(genesDTout) >0){
      message("quired ",length(genes), " genes, finally ",nrow(genesDTout)," records matched!")
      return(genesDTout)
    }else{
      message("quired ",length(genes), " genes, finally 0 records matched!\nplease check your input genes!")
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
    genesDTout <- genesDTout[,-c("geneSymbolUpper")]
    if( nrow(genesDTout) >0){
      message("quired ",length(genes), " genes, finally ",nrow(genesDTout)," records matched!")
      return(genesDTout)
    }else{
      message("quired ",length(genes), " genes, finally 0 records matched!\nplease check your input genes!")
      return(data.table::data.table(genes=genes))
    }
    # entrezId
  }else if( geneType=="entrezId" ){
    if( !all(unlist(lapply(genes, is.numeric))) ){
      stop("Integer is required for entrezId!")
    }
    genesDT <- data.table::data.table(genes=genes)
    gencodeGeneInfo$genes <- gencodeGeneInfo$entrezGeneId
    genesDTout <- merge(genesDT, gencodeGeneInfo, by ="genes", sort=FALSE)
    if( nrow(genesDTout) >0){
      message("quired ",length(genes), " genes, finally ",nrow(genesDTout)," records matched!")
      return(genesDTout)
    }else{
      message("quired ",length(genes), " genes, finally 0 records matched!\nplease check your input genes!")
      return(data.table::data.table(genes=genes))
    }
  }else{
    return(data.table::data.table(genes=genes))
  }
}

#' @title  Fetch information of samples used in analyses from all datasets.
#'
#' @param tissueSiteDetail must be following terms:
#'  \enumerate{
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
#' @return
#' @export
#'
#' @examples
#'  GTExquery_sample( tissueSiteDetail="Liver", dataType="RNASEQ", datasetId="gtex_v8",pageSize=200 )
#'  GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ", datasetId="gtex_v8",pageSize=200 )
#'  GTExquery_sample( tissueSiteDetail="Adipose - Visceral (Omentum)", dataType="RNASEQ", datasetId="gtex_v8",pageSize=200 )
GTExquery_sample <- function( tissueSiteDetail="Liver", dataType="RNASEQ", datasetId="gtex_v8", pageSize=200 ){
  page_tmp <- 0
  pageSize_tmp <- as.integer(pageSize)

  # parameter check:
  if( pageSize_tmp > 2000 | pageSize_tmp < 1){
    stop("pageSize must between 1 and 2000!")
  }
  url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                 "datasetId=", datasetId,"&",
                 "tissueSiteDetail=", tissueSiteDetail,"&",
                 "dataType=",dataType,"&",
                 "page=",page_tmp,"&",
                 "pageSize=", pageSize_tmp, "&",
                 "sortBy=sampleId&sortDirection=asc"
  )
  url1 <- utils::URLencode(url1)
  # url1 <- "https://gtexportal.org/rest/v1/dataset/sample?datasetId=gtex_v8&tissueSiteDetail=Bladder&dataType=RNASEQ&format=json&page=0&pageSize=2000&sortBy=sampleId&sortDirection=asc"
  outInfo <- data.table::data.table()
  pingOut <- apiAdmin_ping()
  if( !is.null(pingOut) & pingOut==200 ){
    message("GTEx API successfully accessed!")
    url1Get <- curl::curl_fetch_memory(url1)
    url1GetText <- rawToChar(url1Get$content)
    url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
    tmp <- data.table::as.data.table(url1GetText2Json$sample)
    outInfo <- rbind(outInfo, tmp)
    message("Total records: ",url1GetText2Json$recordsFiltered,"; total pages: ",url1GetText2Json$numPages,"; downloaded pages:",url1GetText2Json$page," records: ", nrow(tmp))
    while(url1GetText2Json$page < url1GetText2Json$numPages){
      page_tmp <- page_tmp+1
      url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                     "datasetId=", datasetId,"&",
                     "tissueSiteDetail=", tissueSiteDetail,"&",
                     "dataType=",dataType,"&",
                     "page=",page_tmp,"&",
                     "pageSize=", pageSize_tmp, "&",
                     "sortBy=sampleId&sortDirection=asc"
      )
      url1 <- utils::URLencode(url1)
      url1Get <- curl::curl_fetch_memory(url1)
      url1GetText <- rawToChar(url1Get$content)
      url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
      tmp <- data.table::as.data.table(url1GetText2Json$sample)
      outInfo <- rbind(outInfo, tmp)
      message("Total records: ",url1GetText2Json$recordsFiltered,"; total pages: ",url1GetText2Json$numPages,"; downloaded pages:",url1GetText2Json$page," records: ", nrow(tmp))
    }
  }else{
    message("API server is overloaded; please wait several minutes!")
  }
  # note: pathologyNotes info is ignored.
  outInfo[,.(sampleId, sex, ageBracket,datasetId, tissueSiteDetail, tissueSiteDetailId, pathologyNotes, hardyScale,  dataType )]
}

url1 <- "https://gtexportal.org/rest/v1/dataset/sample?datasetId=gtex_v8&tissueSiteDetailId=Liver&dataType=RNASEQ&format=json&page=0&pageSize=2000&sortBy=sampleId&sortDirection=asc"
url1 <- "https://gtexportal.org/rest/v1/dataset/sample?datasetId=gtex_v8&tissueSiteDetail=All&dataType=RNASEQ&format=json&page=0&pageSize=200&sortBy=sampleId&sortDirection=asc"


