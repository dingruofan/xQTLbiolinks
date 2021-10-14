#' @title download gtf files and extract gene info from attribute column:
#'
#' @param gencodeVersion "v26" or "v19"
#'
#' @return
#'
#' @importFrom utils download.file
#' @importFrom data.table fread rbindlist setnames as.data.table data.table
#' @importFrom stringr str_split
#' @importFrom usethis use_data
#' @examples
#' \dontrun{
#'   gtfSubsGeneInfo("v26")
#'   gtfSubsGeneInfo("v19")
#'  }
gtfSubsGeneInfo <- function(gencodeVersion="v26"){
  gtfDir <- tempdir()
  dir.create(gtfDir, recursive = TRUE)
  message("created temp dir: ",gtfDir)
  if(gencodeVersion=="v26" ){
    gtfUrl <- paste0("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz")
    utils::download.file(gtfUrl,paste0(gtfDir,"/gencode.annotation.gtf.gz"), method="curl")
    # downloader::download(gtfUrl,paste0(gtfDir,"/gencode.annotation.gtf.gz"))
  }
  if(gencodeVersion=="v19"){
    gtfUrl <- paste0("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz")
    utils::download.file(gtfUrl,paste0(gtfDir,"/gencode.annotation.gtf.gz"), method="curl")
    # downloader::download(gtfUrl,paste0(gtfDir,"/gencode.annotation.gtf.gz"))
  }
  # deal with gtf:
  if( file.exists(paste0(gtfDir,"/gencode.annotation.gtf.gz")) ){
    message(gtfDir,"/gencode.annotation.gtf.gz", " exists!")
    # read gencode file:
    gencodeAnno <- data.table::fread(paste0(gtfDir,"/gencode.annotation.gtf.gz"),sep="\t", header = FALSE)
    data.table::setnames(gencodeAnno, names(gencodeAnno), c("chr","source","type","start","end","score","strand","phase","attributes") )
    # for gene:
    gencodeAnnoGene <- gencodeAnno[type == "gene"]
# use internal function: gtfSubsGene
    gencodeENSG <- data.table::rbindlist(lapply(gencodeAnnoGene$attributes, gtfSubsGene, att_of_interest= c("gene_id", "gene_type", "gene_name")))
    # gencodeENSG$ENSG <- unlist(lapply(gencodeENSG$gene_id, function(x){ str_split(x, fixed("."))[[1]][1] }))
    if(gencodeVersion=="v26"){
      rowRange1 <- c(seq(from=1, to=nrow(gencodeENSG), by=200), nrow(gencodeENSG))
      gencodeGeneInfoV26 <- data.table::data.table()
      for(i in 1:(length(rowRange1)-1)){
# use internal function: apiRef_gene
        apiRef_geneOut_tmp = apiRef_gene( gencodeENSG[rowRange1[i]:(rowRange1[i+1]-1),]$gene_id, "v26", "GRCh38/hg38")
        gencodeGeneInfoV26 <- rbind(gencodeGeneInfoV26, apiRef_geneOut_tmp)
        rm(apiRef_geneOut_tmp)
        message(i,"/",(length(rowRange1)-1), "; Complete ", round(i/(length(rowRange1)-1)*100,2),"%")
      }
      # add last row:
      apiRef_geneOut_tmp = apiRef_gene( gencodeENSG[rowRange1[length(rowRange1)],]$gene_id, "v26", "GRCh38/hg38")
      gencodeGeneInfoV26 <- rbind(gencodeGeneInfoV26, apiRef_geneOut_tmp)
      rm(apiRef_geneOut_tmp)
      # create .rds file with gencodeGeneInfoV26
      usethis::use_data(gencodeGeneInfoV26)
    }
    if(gencodeVersion=="v19"){
      rowRange1 <- c(seq(from=1, to=nrow(gencodeENSG), by=200), nrow(gencodeENSG))
      gencodeGeneInfoV19 <- data.table::data.table()
      for(i in 1:(length(rowRange1)-1)){
# use internal function: apiRef_gene
        apiRef_geneOut_tmp = apiRef_gene( gencodeENSG[rowRange1[i]:(rowRange1[i+1]-1),]$gene_id, "v19", "GRCh37/hg19")
        gencodeGeneInfoV19 <- rbind(gencodeGeneInfoV19, apiRef_geneOut_tmp)
        rm(apiRef_geneOut_tmp)
        message(i,"/",(length(rowRange1)-1), "; Complete ", round(i/(length(rowRange1)-1)*100,2),"%")
      }
      # add last row:
      apiRef_geneOut_tmp = apiRef_gene( gencodeENSG[rowRange1[length(rowRange1)],]$gene_id, "v19", "GRCh37/hg19")
      gencodeGeneInfoV19 <- rbind(gencodeGeneInfoV19, apiRef_geneOut_tmp)
      rm(apiRef_geneOut_tmp)
      # create .rds file with gencodeGeneInfoV19
      usethis::use_data(gencodeGeneInfoV19, overwrite = TRUE)
    }
  }
}


#' @title Extract gene attributes of interest
#' @description
#' as a funciton of lapply
#' @param gtf_attributes like: c("gene_id", "gene_type", "gene_name")
#'
#' @return specificed attributes
#' @importFrom data.table as.data.table
#' @examples
#'  \dontrun{
#'   # extract gene info:
#'   gencodeENSG <- data.table::rbindlist(lapply(gencodeAnnoGene$attributes, gtfSubsGene, att_of_interest= c("gene_id", "gene_type", "gene_name")))
#'   # extract transcript info:
#'   gencodeENSG <- data.table::rbindlist(lapply(gencodeAnnoGene$attributes, gtfSubsGene, att_of_interest= c("gene_id","transcript_id", "gene_type", "gene_name")))
#'  }
gtfSubsGene <- function(gtf_attributes,  att_of_interest= c("gene_id", "gene_type", "gene_name")){
  att <- unlist(stringr::str_split(gtf_attributes, " ")[[1]])
  # att_of_interest <- c("gene_id", "gene_type", "gene_name")
  if(any(att_of_interest %in% att)){
    dt <- data.table::as.data.table(matrix(gsub("\"|;","", att[which(att %in% att_of_interest)+1]),nrow=1))
    names(dt) <- att_of_interest
    return( dt )
  }else{
    return(NA)}
}


#' @title fetch reference genes by API.
#'
#' @param geneId a character or a character vector(versioned ensemble ID or unversioned ensemble ID)
#' @param gencodeVersion "v26" or "v19"
#' @param genomeBuild "GRCh38/hg38" or "GRCh37/hg19"
#'
#' @return queried gene info table
#' @importFrom jsonlite fromJSON
#' @importFrom data.table as.data.table
#' @export
#' @examples
#' \dontrun{
#'   apiRef_gene(c("ENSG00000116885.18","ENSG00000222623"), "v26", "GRCh38/hg38")
#'  }
apiRef_gene <- function(geneId="", gencodeVersion="v26", genomeBuild="GRCh38/hg38" ){
  geneId <- as.character(unlist(geneId))
  if( any(is.na(geneId)) | any(geneId=="") | any(is.null(geneId)) ){
    stop("gene ID can not be null!")
  }else if(length(geneId)>200){
    stop("number of gene ID can not > 200 ")
  }else{
    url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                   "geneId=", paste0(geneId, collapse = ","),"&",
                   "gencodeVersion=", gencodeVersion,"&",
                   "genomeBuild=",genomeBuild,"&",
                   "page=0&pageSize=2000&format=json"
                   )
# use internal function: apiAdmin_ping
    pingOut <- apiAdmin_ping()
    if( !is.null(pingOut) && pingOut==200 ){
      message("GTEx API successfully accessed!")
      # url1Get <- httr::GET(url1, httr::progress())
      url1Get <- curl::curl_fetch_memory(url1)
      # url1GetText <- httr::content(url1Get,"text", encoding = "UTF-8")
      url1GetText <- rawToChar(url1Get$content)
      url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
      url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
      url1GetText2Json2DT$genomeBuild <- genomeBuild
      outInfo <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
      return(outInfo)
    }else{
      stop("GTEx API can not be accessed, please check your network!")
    }
  }
}


#' @title create .rds file with GTExquery_sample function
#'
#' @param datasetId "gtex_v8" or "gtex_v7"
#' @import usethis
#'
#' @return none
#'
#' @examples
#' \dontrun{
#'   createTissueSiteDetailMappingData("gtex_v8")
#'   createTissueSiteDetailMappingData("gtex_v7")
#'  }
createTissueSiteDetailMappingData <- function(datasetId="gtex_v8"){
  # obtain all tissueSiteDetail info:
  if( datasetId == "gtex_v8" ){
    tissueSiteDetailGTExv8 <- GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ", datasetId="gtex_v8",pageSize=2000 )
    tissueSiteDetailGTExv8 <- unique(tissueSiteDetailGTExv8[,.(tissueSiteDetail,tissueSiteDetailId)][order(tissueSiteDetail)])
    usethis::use_data(tissueSiteDetailGTExv8, overwrite = TRUE)
  }else if(datasetId == "gtex_v7" ){
    tissueSiteDetailGTExv7 <- GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ", datasetId="gtex_v7",pageSize=2000 )
    tissueSiteDetailGTExv7 <- unique(tissueSiteDetailGTExv7[,.(tissueSiteDetail,tissueSiteDetailId)][order(tissueSiteDetail)])
    usethis::use_data(tissueSiteDetailGTExv7, overwrite = TRUE)
  }
}


# httr::use_proxy(url="127.0.0.1", port=7890
#           # ,username="dd",password="123456"
# )

# usethis::use_package("data.table")
# usethis::use_package("curl")
# usethis::use_package("jsonlite")
# usethis::use_package("stringr")
# usethis::use_package("usethis")
# usethis::use_package("utils")






