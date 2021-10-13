#' Title download gtf files and extract gene info from attribute column:
#'
#' @param gencodeVersion
#'
#' @return
#'
#' @importFrom utils download.file
#' @importFrom data.table fread rbindlist setnames as.data.table
#' @importFrom stringr str_split
#' @examples
#'  gtfSubsGeneInfo("v26")
gtfSubsGeneInfo <- function(gencodeVersion="v26"){
  gtfDir <- tempdir()
  dir.create(gtfDir, recursive = TRUE)
  message("created temp dir: ",gtfDir)
  if(gencodeVersion=="v26"){
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
# use internal function: apiRef_gene
    if(gencodeVersion=="v26"){
      rowRange1 <- c(seq(from=1, to=nrow(gencodeENSG), by=1000), nrow(gencodeENSG))
      for(i in 1:(length(rowRange1)-1)){
        a = apiRef_gene(gencodeENSG[rowRange1[i]:(rowRange1[i+1]-1),]$gene_id, "v26", "GRCh38/hg38")
      }
    }
  }
}


#' Title extract gene attributes of interest
#' @description
#' as a funciton of lapply
#' @param gtf_attributes
#'
#' @return specificed attributes
#' @examples
#'  # extract gene info:
#'  gencodeENSG <- data.table::rbindlist(lapply(gencodeAnnoGene$attributes, gtfSubsGene, att_of_interest= c("gene_id", "gene_type", "gene_name")))
#'  # extract transcript info:
#'  gencodeENSG <- data.table::rbindlist(lapply(gencodeAnnoGene$attributes, gtfSubsGene, att_of_interest= c("gene_id","transcript_id", "gene_type", "gene_name")))
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


#' Title fetch reference genes by API.
#'
#' @param geneId a character or a character vector(versioned ensemble ID or unversioned ensemble ID)
#' @param gencodeVersion "v26" or "v19"
#' @param genomeBuild "GRCh38/hg38" or "GRCh37/hg19"
#'
#' @return queried gene info table
#' @importFrom httr GET progress content
#' @importFrom jsonlite fromJSON
#' @importFrom data.table as.data.table
#' @examples
#'  apiRef_gene(c("ENSG00000116885.18","ENSG00000222623"), "v26", "GRCh38/hg38")
apiRef_gene <- function(geneId="", gencodeVersion="v26", genomeBuild="GRCh38/hg38" ){
  geneId <- as.character(unlist(geneId))
  if( any(is.na(geneId)) | any(geneId=="") | any(is.null(geneId)) ){
    stop("gene ID can not be null!")
  }else if(length(geneId)>2000){
    stop("number of gene ID can not > 2000 ")
  }else{
    url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                   "geneId=", paste0(geneId, collapse = ","),"&",
                   "gencodeVersion=", gencodeVersion,"&",
                   "genomeBuild=",genomeBuild,"&",
                   "page=0&pageSize=2000&format=json"
                   )
# use internal function: apiAdmin_ping
    if( apiAdmin_ping()==200 ){
      message("GTEx API successfully accessed!")
      url1Get <- httr::GET(url1, httr::progress())
      url1GetText <- httr::content(url1Get,"text", encoding = "UTF-8")
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

#
#' @title Heartbeat to check server connectivity.
#' @description
#'  test API server
#' @importFrom httr GET status_code
#' @return boolean value
#' @examples
#'  apiAdmin_ping()
#'  apiStatus <- ifelse( apiAdmin_ping() ==200, "GTEx API can be accessed", "Please check your network!")
#'  print(apiStatus)
apiAdmin_ping <- function(){
  url1Get <- "https://gtexportal.org/rest/v1/admin/ping"
  tryCatch(
    {
      httr::status_code(httr::GET(url1Get))
    },
    # e = simpleError("test error"),
    error=function(cond){
      message(cond,"\nplease check your network!")
    },
    warning = function(cond){
      message(cond)
    }
  )
}

# httr::use_proxy(url="127.0.0.1", port=7890
#           # ,username="dd",password="123456"
# )

