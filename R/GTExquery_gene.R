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

#' @Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
GTExquery_sample <- function( tissueSiteDetailId="Liver", dataType="RNASEQ", datasetId="gtex_v8" ){
  page_tmp <- 0
  pageSize_tmp <- 200
  url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                 "datasetId=", datasetId,"&",
                 "tissueSiteDetailId=", tissueSiteDetailId,"&",
                 "dataType=",dataType,"&",
                 "page=",page_tmp,"&",
                 "pageSize=", pageSize_tmp, "&",
                 "sortBy=sampleId&sortDirection=asc"
  )
  # url1 <- "https://gtexportal.org/rest/v1/dataset/sample?datasetId=gtex_v8&tissueSiteDetailId=Bladder&dataType=RNASEQ&format=json&page=0&pageSize=2000&sortBy=sampleId&sortDirection=asc"
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
                     "tissueSiteDetailId=", tissueSiteDetailId,"&",
                     "dataType=",dataType,"&",
                     "page=",page_tmp,"&",
                     "pageSize=", pageSize_tmp, "&",
                     "sortBy=sampleId&sortDirection=asc"
      )
      url1Get <- curl::curl_fetch_memory(url1)
      url1GetText <- rawToChar(url1Get$content)
      url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
      tmp <- data.table::as.data.table(url1GetText2Json$sample)
      outInfo <- rbind(outInfo, tmp)
      message("Total records: ",url1GetText2Json$recordsFiltered,"; total pages: ",url1GetText2Json$numPages,"; downloaded pages:",url1GetText2Json$page," records: ", nrow(tmp))
    }
  }
  # note: pathologyNotes info is ignored.
  outInfo[,.(sampleId, sex, ageBracket,datasetId, tissueSiteDetail, tissueSiteDetailId, pathologyNotes, hardyScale,  dataType )]
}

url1 <- "https://gtexportal.org/rest/v1/dataset/sample?datasetId=gtex_v8&tissueSiteDetailId=Bladder&format=json&page=0&pageSize=200&sortBy=sampleId&sortDirection=asc"


