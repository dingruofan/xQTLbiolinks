#' @title Query gene information through API.
#' @description
#'  users can use this function to search gene of interest
#' @param genes A gene symbol, gencode id (versioned), or a charater string of gene type.
#' \itemize{
#'   \item \strong{gene symbol (Default)}.
#'
#'   A character string or a character vector (case ignored). like: "tp53","naDK","SDF4".
#'   \item \strong{gencode/ensemble id} (versioned or unversioned).
#'
#'    A character string or a character vector (case ignored). like: "ENSG00000210195.2","ENSG00000078808".

#'   \item \strong{gene classification}.
#'
#'   "genes" must be chosen from following list when "geneType" is "geneCategory".
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
#' @param geneType A character string in "geneSymbol"(default), "gencodeId" and "geneCategory".
#'
#' @param gencodeVersion "v26" or "v19"
#' @param recordPerChunk A integer value. Defaulut: 150
#'
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
#' @import data.table
#' @import stringr
#' @import jsonlite
#' @import utils
#' @import curl
#' @import stats
#' @export
#'
#' @examples
#' \donttest{
#'  # get all protein coding genes and description:
#'  protein_coding <- GTExquery_gene(genes="protein coding", geneType="geneCategory", "v26" )
#'  # get all miRNA and description:
#'  miRNA <- GTExquery_gene(genes="miRNA", geneType="geneCategory", "v19")
#'
#'  # hg38 test:
#'  geneInfo <- GTExquery_gene("TP53", "geneSymbol", "v26")
#'  geneInfo <- GTExquery_gene(c("tp53","naDK","SDF4"), "geneSymbol", "v26")
#'  geneInfo <- GTExquery_gene(c("ENSG00000210195.2","ENSG00000078808"),
#'                               geneType="gencodeId", "v26")
#'  # hg19 test:
#'  geneInfo <- GTExquery_gene(c("TP53","naDK"), "geneSymbol", "v19")
#'  geneInfo <- GTExquery_gene(c("ENSG00000141510.11","ENSG00000008130.11"), "gencodeId", "v19")
#'  }
GTExquery_gene <- function(genes="", geneType="geneSymbol", gencodeVersion="v26", recordPerChunk=150){
  geneSymbol <- gencodeId <- entrezGeneId <- chromosome <- start <- end <- strand <- tss <- description <- cutF <- genesUpper <- NULL
  .<-NULL
  page_tmp <- 0
  pageSize_tmp <- recordPerChunk
  cutNum <- recordPerChunk
  genomeBuild="GRCh38/hg38"

  # check genes
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }

  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("geneSymbol", "gencodeId", "geneCategory")) ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }

  # gencodeVersion
  if( is.null(gencodeVersion) ||  any(is.na(gencodeVersion)) || any(gencodeVersion=="") || length(gencodeVersion)!=1){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }else if( !(gencodeVersion %in% c("v26", "v19")) ){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }
  # set genomeBuild:
  if(gencodeVersion == "v26"){
    genomeBuild="GRCh38/hg38"
  }else if(gencodeVersion == "v19"){
    genomeBuild="GRCh37/hg19"
  }

  ######################### if geneType is "geneCategory":
  if( geneType == "geneCategory" ){
    # Fetch all genes' info:
    gencodeGeneInfo <- GTExquery_geneAll("v26")
    geneCategory = unique(gencodeGeneInfo$geneType)
    if( length(genes)!=1 || any(!(genes %in% geneCategory)) ){
      message(paste0(1:length(geneCategory),". ",unique(rev(gencodeGeneInfo$geneType)), collapse = "\n"))
      stop("if \"geneType\" is  \"geneCategory\", input \"genes\" should be choosen from the above: ")
    }else{
      genesDTout <- cbind(data.table::data.table(genes=gencodeGeneInfo[geneType==genes,]$geneSymbol), gencodeGeneInfo[geneType==genes,])
      return(genesDTout)
    }
  }else{
    ######################### if geneType is "geneSymbol" or "genecodeId":
    # API处理GENE，单个字符串模糊匹配，多个字符串精确匹配:
    if( length(genes)==1 ){
      # construct url:
      if( stringr::str_detect(genes,stringr::regex("^ENSG00000")) ){
        url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                       "geneId=", genes,"&",
                       "gencodeVersion=", gencodeVersion,"&",
                       "genomeBuild=",genomeBuild,"&",
                       "page=",page_tmp,"&",
                       "pageSize=", pageSize_tmp,"&",
                       "format=json"
        )
      }else{
        url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                       "geneId=^", genes,"$&",
                       "gencodeVersion=", gencodeVersion,"&",
                       "genomeBuild=",genomeBuild,"&",
                       "page=",page_tmp,"&",
                       "pageSize=", pageSize_tmp,"&",
                       "format=json"
        )
      }
      url1 <- utils::URLencode(url1)
      outInfo <- data.table::data.table()
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
        if( nrow(url1GetText2Json2DT)==0 ){
          message( "0 record fatched!" )
          return(data.table::data.table())
        }else{
          url1GetText2Json2DT$genomeBuild <- genomeBuild
          tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
          outInfo <- rbind(outInfo, tmp)
          outInfo <- cbind(data.table(genes=genes), outInfo)
          message( nrow(tmp), " record has been obtained!" )
          return(outInfo)
        }
      }else{
        message("")
        return(data.table::data.table())
      }
    }else if( length(genes)>1 ){
      # 分批下载：
      genesCut <- data.table::data.table(genes=genes, ID=1:length(genes), cutF = as.character(cut(1:length(genes),breaks=seq(0,length(genes)+cutNum,cutNum) )) )
      genesCut$genesUpper <- toupper(genesCut$genes)
      genesURL <- genesCut[,.(genesURL=paste0(genes,collapse = "%2C")),by=c("cutF")]
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
        url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                       "geneId=", genesURL[i,]$genesURL,"&",
                       "gencodeVersion=", gencodeVersion,"&",
                       "genomeBuild=",genomeBuild,"&",
                       "page=",page_tmp,"&",
                       "pageSize=", pageSize_tmp,"&",
                       "format=json"
        )
        url1 <- utils::URLencode(url1)
        url1Get <- curl::curl_fetch_memory(url1)
        url1GetText <- rawToChar(url1Get$content)
        url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
        url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
        if( nrow(url1GetText2Json2DT)==0 ){
          message( "0 record fatched!" )
          return(data.table::data.table())
        }
        url1GetText2Json2DT$genomeBuild <- genomeBuild
        tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
        tmp$genesUpper <- toupper(unlist(tmp[,geneType,with=FALSE]))
        # because of versioned and unversioned gencodeID, merge separately is needed!
        if(geneType == "gencodeId"){
          # versioned:
          tmp1 <- merge( genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp, by="genesUpper", sort = FALSE)[,-c("genesUpper")]
          # unversioned:
          tmp$genesUpper <- unlist(lapply(tmp$genesUpper, function(x){ stringr::str_split(x,stringr::fixed("."))[[1]][1] }))
          tmp2 <- merge( genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp, by="genesUpper", sort = FALSE)[,-c("genesUpper")]
          # combination:
          tmp <- merge(genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], rbind(tmp1, tmp2), by="genes",all.x=TRUE,sort = FALSE)
        }else{
          tmp <- merge( genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp, by="genesUpper",all.x=TRUE, sort = FALSE)[,-c("genesUpper")]
        }
        outInfo <- rbind(outInfo, tmp)
        message("Downloaded  ", i, "/",nrow(genesURL),"; ", length(na.omit(outInfo$gencodeId)), " records.")
        # message("Downloaded  ", round(i/nrow(genesURL)*100,2),"%; totally ", length(na.omit(outInfo$gencodeId)), " records fetched!")
        # rm(url1, url1Get, url1GetText, url1GetText2Json, url1GetText2Json2DT)
      }
      return(outInfo)
    }else{
      return(data.table::data.table())
    }
  }
}

#' @title  Fetch information of samples used in analyses from all datasets.
#'
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
#' @param dataType A character string. Options: "RNASEQ" (default), "WGS", "WES", "OMNI".
#' @param datasetId A character string. Options: "gtex_v8" (default), "gtex_v7".
#' @param recordPerChunk A integer value (1-2000). number of records fetched per request (default: 200).
#' @import data.table
#' @import utils
#' @import curl
#' @import jsonlite
#' @return returen sample information
#' @export
#' @examples
#' \donttest{
#'   GTExquery_sample( tissueSiteDetail="Liver", dataType="RNASEQ",
#'                     datasetId="gtex_v8", recordPerChunk=200 )
#'   GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ",
#'                     datasetId="gtex_v8",recordPerChunk=2000 )
#'   GTExquery_sample( "Adipose - Visceral (Omentum)", "RNASEQ", "gtex_v8", 200 )
#'   }
GTExquery_sample <- function( tissueSiteDetail="Liver", dataType="RNASEQ", datasetId="gtex_v8", recordPerChunk=200 ){
  sampleId <- sex <- ageBracket <- pathologyNotes <- hardyScale <- NULL
  .<-NULL
  page_tmp <- 0
  pageSize_tmp <- as.integer(recordPerChunk)

  ########## parameter check: tissueSiteDetail
  if( is.null(datasetId) ||  any(is.na(datasetId)) ){
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
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail))  ){
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
  if(is.null(pageSize_tmp) ||  any(is.na(pageSize_tmp)) ){
    stop("Parameter \"pageSize\" should be a integer!")
    }else if(length(pageSize_tmp)!=1){
      stop("Parameter \"pageSize\" should be a integer!")
    }else if( pageSize_tmp > 2000 | pageSize_tmp < 1){
      stop("pageSize must between 1 and 2000!")
    }

  ######## parameter check: dataType
  if(is.null(dataType) ||  any(is.na(dataType))){
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
    message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp<-page_tmp+1
    while( page_tmp <= (url1GetText2Json$numPages-1)){
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
      message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
      page_tmp <- page_tmp+1
    }
    return(outInfo[,.(sampleId, sex, ageBracket,datasetId, tissueSiteDetail, tissueSiteDetailId, pathologyNotes, hardyScale,  dataType )])

  }else{
    message("")
  }
  # note: pathologyNotes info is ignored.
}

# Indepedently fetch geneInfo:
#' Title
#'
#' @param gencodeVersion A character string. "v26" or "v19".
#' @param recordPerChunk A integer value. From 1 to 2000.
#'
#' @return A data.table with all genes' information.
#' @import utils
#' @import data.table
GTExquery_geneAll <- function(gencodeVersion="v26", recordPerChunk=2000){
  geneSymbol <- gencodeId <- entrezGeneId <- geneType <- chromosome <- start <- end <- strand <- tss <- description <- NULL
  .<-NULL
  page_tmp <- 0
  pageSize_tmp <- recordPerChunk
  genomeBuild="GRCh38/hg38"
  # gencodeVersion
  if( is.null(gencodeVersion) ||  any(is.na(gencodeVersion)) || any(gencodeVersion=="") || length(gencodeVersion)!=1){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }else if( !(gencodeVersion %in% c("v26", "v19")) ){
    stop("Parameter \"gencodeVersion\" should be choosen from \"v26\", \"v19\".")
  }
  # set genomeBuild:
  if(gencodeVersion == "v26"){
    genomeBuild="GRCh38/hg38"
  }else if(gencodeVersion == "v19"){
    genomeBuild="GRCh37/hg19"
  }

  # construct url:
  url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                 "geneId=", "[a-zA-Z0-9]","&",
                 "gencodeVersion=", gencodeVersion,"&",
                 "genomeBuild=",genomeBuild,"&",
                 "page=",page_tmp,"&",
                 "pageSize=", pageSize_tmp,"&",
                 "format=json"
  )
  url1 <- utils::URLencode(url1)
  outInfo <- data.table::data.table()
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
    tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
    outInfo <- rbind(outInfo, tmp)
    message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp <- page_tmp+1
    while( page_tmp <= (url1GetText2Json$numPages-1) ){
      url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                     "geneId=", "[a-zA-Z0-9]","&",
                     "gencodeVersion=", gencodeVersion,"&",
                     "genomeBuild=",genomeBuild,"&",
                     "page=",page_tmp,"&",
                     "pageSize=", pageSize_tmp,"&",
                     "format=json"
      )
      url1 <- utils::URLencode(url1)
      url1Get <- curl::curl_fetch_memory(url1)
      url1GetText <- rawToChar(url1Get$content)
      url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
      url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
      url1GetText2Json2DT$genomeBuild <- genomeBuild
      tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
      outInfo <- rbind(outInfo, tmp)
      message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
      page_tmp <- page_tmp+1
    }
    return(outInfo)
  }else{
    stop("Network error!")
    return(data.table::data.table())
  }
}


#' @title Query variant information by snpId or variant ID.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param variantType A character string. "snpId" or "variantId". Default: "snpId".
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#'
#' @return A data.table.
#' @export
#'
#' @examples
#'  \donttest{
#'   GTExquery_varId("rs12596338")
#'   GTExquery_varId("rs12596338", datasetId="gtex_v7")
#'   GTExquery_varId("chr11_66561248_T_C_b38", variantType="variantId", datasetId="gtex_v8")
#'   GTExquery_varId("11_66328719_T_C_b37", variantType="variantId", datasetId="gtex_v7")
#'  }
GTExquery_varId <- function(variantName="", variantType="snpId", datasetId="gtex_v8"){
  ########## parameter check: variantName
  if(is.null(variantName) ||  any(is.na(variantName)) ){
    stop("Parameter \"variantName\" can not be NULL or NA!")
  }else if(length(variantName)!=1){
    stop("Parameter \"variantName\" can not be NULL or NA!")
  }else if( variantType=="snpId" && !(stringr::str_detect(variantName, stringr::regex("^rs")))  ){
    stop("Parameter \"variantName\" must begin with a \"rs\", like: \"rs147502335\", \"rs147538909\".")
  }
  ########## parameter check: variantType
  if(is.null(variantType) ||  any(is.na(variantType)) ){
    stop("Parameter \"variantType\" can not be NULL or NA!")
  }else if(length(variantType)!=1){
    stop("Parameter \"variantType\" can not be NULL or NA!")
  }else if(  !(variantType %in% c("snpId", "variantId"))  ){
    stop("Parameter \"variantType\" should be chosen from \"snpId\" or \"variantId\"!")
  }
  ########## parameter check: datasetId
  if(is.null(datasetId) ||  any(is.na(datasetId)) ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(  !(datasetId %in% c("gtex_v8", "gtex_v7"))  ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }

  ############# if variantType is "snpId":
  if( variantType == "snpId" ){
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/variant??format=json","&",
                   "datasetId=", datasetId,"&",
                   "snpId=", variantName)
    ############# if variantType is "variantId":
  }else if( variantType == "variantId" ){
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/variant??format=json","&",
                   "datasetId=", datasetId,"&",
                   "variantId=", variantName)
  }
  # url1 <- "https://gtexportal.org/rest/v1/dataset/variant?format=json&datasetId=gtex_v8&snpId=rs12596338"
  url1 <- utils::URLencode(url1)
  # test api server accessibility:
  pingOut <- apiAdmin_ping()
  if( !(!is.null(pingOut) && pingOut==200) ){
    return(data.table::data.table())
  }
  # fetch data:
  url1Get <- curl::curl_fetch_memory(url1)
  url1GetText <- rawToChar(url1Get$content)
  url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
  tmp <- data.table::as.data.table(url1GetText2Json$variant)
  if(nrow(tmp)==0){
    message("No variant found, please check your input.")
    return(data.table::data.table())
  }
  # eliminate the difference between the gtex_v7 and gtex_v8:
  if(datasetId == "gtex_v7"){
    tmp$b37VariantId <- tmp$variantId
  }
  outInfo <- tmp[,c("variantId", "snpId","b37VariantId", "chromosome","pos", "ref", "alt","datasetId","maf01", "shorthand")]
  return(outInfo)
}

#' @title Query variant information by position.
#'
#' @param chrom A character string.
#' @param pos An integer array.
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#'
#' @return A data.table.
#' @export
#'
#' @examples
#' \donttest{
#'  GTExquery_varPos(chrom="chr1", pos=c(1102708,1105739),"gtex_v8")
#'  GTExquery_varPos(chrom="1", pos=c(1038088,1041119),"gtex_v7")
#' }
GTExquery_varPos <- function(chrom="", pos=numeric(0), datasetId="gtex_v8"){
  pos
  ########## parameter check: variantName
  if(is.null(chrom) ||  any(is.na(chrom)) ){
    stop("Parameter \"chrom\" can not be NULL or NA!")
  }else if(length(chrom)!=1){
    stop("Parameter \"chrom\" can not be NULL or NA!")
  }else if( datasetId=="gtex_v8" && !(chrom %in% c(paste0("chr",c(1:23,"X","Y"))))  ){
    stop("For GTEx v8, Parameter \"chrom\" must begin with a \"chr\", like: \"chr1\", \"chrX\".")
  }else if( datasetId=="gtex_v7" && !(chrom %in% c(1:23,"X","Y"))  ){
    stop("For GTEx v7, Parameter \"chrom\" must begin without a \"chr\", like: \"1\", \"X\".")
  }
  ########## parameter check: pos
  if(is.null(pos) ||  any(is.na(pos)) || length(pos)==0){
    stop("Parameter \"pos\" can not be NULL or NA!")
  }else if( any(!is.wholenumber(pos)) || any(pos<=0) ){
    stop("Parameter \"pos\" should be an positive integer!")
  }
  ########## parameter check: datasetId
  if(is.null(datasetId) ||  any(is.na(datasetId)) ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(  !(datasetId %in% c("gtex_v8", "gtex_v7"))  ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }

  url1 <- paste0("https://gtexportal.org/rest/v1/dataset/variant??format=json","&",
                 "datasetId=", datasetId,"&",
                 "chromosome=", chrom, "&",
                 "pos=",paste0(pos,collapse=","))
  # https://gtexportal.org/rest/v1/dataset/variant?format=json&datasetId=gtex_v8&chromosome=chr11&pos=65592772
  # https://gtexportal.org/rest/v1/dataset/variant?format=json&datasetId=gtex_v7&chromosome=16&pos=57190138
  # https://gtexportal.org/rest/v1/dataset/variant?format=json&datasetId=gtex_v7&chromosome=1&pos=115746%2C135203%2C1086820
  # fetch data:
  url1Get <- curl::curl_fetch_memory(url1)
  url1GetText <- rawToChar(url1Get$content)
  url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
  tmp <- data.table::as.data.table(url1GetText2Json$variant)
  if(nrow(tmp)==0){
    message("No variant found, please check your input.")
    return(data.table::data.table())
  }
  # eliminate the difference between the gtex_v7 and gtex_v8:
  if(datasetId == "gtex_v7"){
    tmp$b37VariantId <- tmp$variantId
  }
  outInfo <- tmp[,c("variantId", "snpId","b37VariantId", "chromosome","pos", "ref", "alt","datasetId","maf01", "shorthand")]
  return(outInfo)
}



#
#' @title Heartbeat to check server connectivity.
#' @description
#'  test API server
#' @import curl
#' @return A numeric value. Response code 200 indicates that the request has succeeded.
#' @export
#' @examples
#' \donttest{
#'   apiAdmin_ping()
#'   ifelse( apiAdmin_ping() ==200, "accessed", "failed")
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


#' @title determine is a whole number:
#'
#' @param x A number
#' @param tol Don't change
#'
#' @return TRUE or FALSE
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}





