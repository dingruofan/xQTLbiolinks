#' @title download gtf files and extract gene info from attribute column:
#'
#' @param gencodeVersion "v26" or "v19"
#'
#' @return create .rds file
#'
#' @importFrom utils download.file
#' @importFrom data.table fread rbindlist setnames as.data.table data.table
#' @importFrom stringr str_split
gtfSubsGeneInfo <- function(gencodeVersion="v26"){
  type <- NULL
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
      # usethis::use_data(gencodeGeneInfoV26)
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
      # usethis::use_data(gencodeGeneInfoV19, overwrite = TRUE)
    }
  }
  return(NULL)
}


#' @title Extract gene attributes of interest
#' @description
#' as a funciton of lapply
#' @param gtf_attributes A character string or a character vector. Like: c("gene_id", "gene_type", "gene_name"). Default: "gene_id", "gene_type", "gene_name".
#' @param att_of_interest A character string or a character vector. Attributes of interest.
#'
#' @return specificed attributes
#' @importFrom data.table as.data.table
#' @examples
#'  \dontrun{
#'   # extract gene info:
#'   gencodeENSG <- data.table::rbindlist(lapply(gencodeAnnoGene$attributes,
#'                                               gtfSubsGene,
#'                                               c("gene_id", "gene_type", "gene_name")))
#'   # extract transcript info:
#'   gencodeENSG <- data.table::rbindlist(lapply(gencodeAnnoGene$attributes,
#'                                               gtfSubsGene,
#'                                               c("gene_id","transcript_id",
#'                                               "gene_type", "gene_name")))
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
  geneSymbol <- gencodeId <- entrezGeneId <- geneType <- chromosome <- start <- end <- strand <- tss <- description <- NULL
  .<-NULL
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
    # check network:
    bestFetchMethod <- apiAdmin_ping()
    if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
      message("Note: API server is busy or your network has latency, please try again later.")
      return(NULL)
    }
    message("GTEx API successfully accessed!")
    # url1Get <- httr::GET(url1, httr::progress())
    url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
    url1GetText2Json2DT$genomeBuild <- genomeBuild
    outInfo <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
    return(outInfo)
  }
}


#' @title create .rds file with GTExquery_sample function
#'
#' @param datasetId "gtex_v8" or "gtex_v7"
#'
#' @return none
createTissueSiteDetailMappingData <- function(datasetId="gtex_v8"){
  # tissueSiteDetail <- tissueSiteDetailId <- NULL
  # .<-NULL
  # # obtain all tissueSiteDetail info:
  # if( datasetId == "gtex_v8" ){
  #   tissueSiteDetailGTExv8 <- GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ", datasetId="gtex_v8",recordPerChunk =2000 )
  #   tissueSiteDetailGTExv8 <- unique(tissueSiteDetailGTExv8[,.(tissueSiteDetail,tissueSiteDetailId)][order(tissueSiteDetail)])
  #   usethis::use_data(tissueSiteDetailGTExv8, overwrite = TRUE)
  # }else if(datasetId == "gtex_v7" ){
  #   tissueSiteDetailGTExv7 <- GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ", datasetId="gtex_v7",recordPerChunk =2000 )
  #   tissueSiteDetailGTExv7 <- unique(tissueSiteDetailGTExv7[,.(tissueSiteDetail,tissueSiteDetailId)][order(tissueSiteDetail)])
  #   usethis::use_data(tissueSiteDetailGTExv7, overwrite = TRUE)
  # }
  return(NULL)
}

#' @title Query gene information through all genes' information
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
#'
#'   \item \strong{geneCategory}.
#'
#'   When choose "geneCategory", "genes" must be chosen from following gene category:
#'   \itemize{
#'   \item protein coding
#'   \item antisense
#'   \item lincRNA
#'   \item unprocessed pseudogene
#'   \item miRNA
#'   \item transcribed unprocessed pseudogene
#'   \item snRNA
#'   \item processed pseudogene
#'   \item processed transcript
#'   \item TEC
#'   \item transcribed unitary pseudogene
#'   \item transcribed processed pseudogene
#'   \item sense intronic
#'   \item misc RNA
#'   \item snoRNA
#'   \item scaRNA
#'   \item rRNA
#'   \item unitary pseudogene
#'   \item 3prime overlapping ncRNA
#'   \item polymorphic pseudogene
#'   \item bidirectional promoter lncRNA
#'   \item sense overlapping
#'   \item pseudogene
#'   \item IG V pseudogene
#'   \item scRNA
#'   \item IG C gene
#'   \item IG J gene
#'   \item IG V gene
#'   \item sRNA
#'   \item ribozyme
#'   \item vaultRNA
#'   \item non coding
#'   \item TR J gene
#'   \item TR C gene
#'   \item TR V gene
#'   \item TR V pseudogene
#'   \item TR D gene
#'   \item IG C pseudogene
#'   \item macro lncRNA
#'   \item TR J pseudogene
#'   \item IG D gene
#'   \item IG J pseudogene
#'   \item IG pseudogene
#'   \item Mt tRNA
#'   \item Mt rRNA
#'   }
#' }
#'
#' @param geneType A character string. Types of queried genes. Options: "geneSymbol" (default), "gencodeId", "entrezId";
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
apiRef_genes <- function(genes="", geneType="geneSymbol", gencodeVersion="v26", genomeBuild="GRCh38/hg38"){
  # check null/na
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }
  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\", \"entrezId\".")
  }else if( !(geneType %in% c("geneSymbol", "gencodeId", "entrezId","geneCategory")) ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\", \"entrezId\",\"geneCategory\".")
  }
  # gencodeVersion
  if( is.null(gencodeVersion) ||  any(is.na(gencodeVersion)) || any(gencodeVersion=="") || length(gencodeVersion)!=1){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }else if( !(gencodeVersion %in% c("v26", "v19")) ){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }
  # genomeBuild
  # check gencodeVersion match with genomeBuild
  if( is.null(genomeBuild) ||  any(is.na(genomeBuild)) || any(genomeBuild=="") || length(genomeBuild)!=1){
    stop("Parameter \"genomeBuild\" should be choosen from \"GRCh38/hg38\", \"GRCh37/hg19\".")
  }else if( !(genomeBuild %in% c("GRCh38/hg38", "GRCh37/hg19")) ){
    stop("Parameter \"genomeBuild\" should be choosen from \"GRCh38/hg38\", \"GRCh37/hg19\".")
  } else if(gencodeVersion=="v26" & genomeBuild=="GRCh38/hg38"){
    # data(gencodeGeneInfoV26)
    # gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV26)
    gencodeGeneInfo <- GTExquery_geneAll("v26")
  }else if(gencodeVersion=="v19" & genomeBuild=="GRCh37/hg19"){
    # data(gencodeGeneInfoV19)
    # gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV19)
    gencodeGeneInfo <- GTExquery_geneAll("v19")
  }else{
    stop("gencodeVersion must be matched with genomeBuild.\n eg. v26(GRCh38/hg38), v19(GRCh37/hg19)")
  }

  # merge:
  if( geneType=="geneSymbol" ){
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
    # geneCategory
  }else if(geneType == "geneCategory"){
    geneCategory = unique(gencodeGeneInfo$geneType)
    if( length(genes)!=1 || any(!(genes %in% geneCategory)) ){
      message(paste0(1:length(geneCategory),". ",geneCategory, collapse = "\n"))
      stop("if \"geneType\" is  \"geneCategory\", input \"genes\" should be choosen from the above: ")
    }else{
      genesDTout <- cbind(data.table(genes=gencodeGeneInfo[geneType==genes,]$geneSymbol), gencodeGeneInfo[geneType==genes,])
      return(genesDTout)
    }
  }else{
    return(data.table::data.table(genes=genes))
  }
}



# httr::use_proxy(url="127.0.0.1", port=7890
#           # ,username="dd",password="123456"
# )


# devtools::build("D:/R_project/GTExbiolinks","D:/R_project/GTExbiolinks.tar.gz")
#
# install.packages("~/GTExbiolinks.tar.gz", repos = NULL, type = "source")
# install.packages("~/GTExbiolinks.tar.gz", repos = NULL, type = "source", lib="/home/dingruofan/anaconda3/envs/work/lib/R/library")
# install.packages("~/GTExbiolinks_0.0.0.9000.tar.gz", repos = NULL, type = "source", lib="/home/dingruofan/anaconda3/envs/work/lib/R/library")

# import:
# usethis::use_package("data.table",min_version ="1.14.2")
# usethis::use_package("curl", min_version = "4.3.2")
# usethis::use_package("jsonlite", min_version = "1.7.2")
# usethis::use_package("stringr", min_version = "1.4.0")
# usethis::use_package("utils", min_version = "4.0.3")
# usethis::use_package("SummarizedExperiment", min_version = "1.18.2")
# usethis::use_package("GenomicRanges", min_version = "1.40.0")
# usethis::use_package("IRanges", min_version = "2.22.2")
# usethis::use_package("GenomeInfoDb", min_version = "1.24.2")
# usethis::use_package("ggplot2", min_version = "3.3.5")
# usethis::use_package("rvest", min_version = "1.0.1")
# usethis::use_package("gridExtra", min_version = "2.3")
# usethis::use_package("tidyr", min_version = "1.1.4")
# usethis::use_package("ggrepel", min_version = "0.9.1")
# usethis::use_package("crayon", min_version = "1.4.1")
# usethis::use_package("httr", min_version = "1.4.2")
# usethis::use_package("RMySQL", min_version = "0.10.22")
# usethis::use_package("DBI", min_version = "1.1.1")

# suggest:
# usethis::use_package("coloc", min_version = "5.1.0", type="Suggests")
# usethis::use_package("usethis", min_version = "2.0.1")
# usethis::use_package("rlang", min_version = "0.4.11")


