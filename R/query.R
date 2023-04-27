#' @title Query basic information (including name, symbol, position and description, etc. ) of genes.
#' @param genes A charater vector or a string of gene symbol, gencode id (versioned or unversioned), or a charater string of gene type.
#' \itemize{
#'   \item \strong{gene symbol (Default)}.
#'
#'   A character string or a character vector (case ignored). like: "tp53","naDK","SDF4".
#'   \item \strong{gencode/ensemble id} (versioned or unversioned).
#'
#'    A character string or a character vector (case ignored). like: "ENSG00000210195.2","ENSG00000078808".

#'   \item \strong{gene classification}.
#'
#'   when "geneType" is "geneCategory", supported "genes" can be listed using function `gencodeGenetype$V26` or `gencodeGenetype$V19`
#'
#' }
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#'
#' @param recordPerChunk (integer) number of records fetched per request (default: 150).
#'
#' @return A data.table object of queried gene information. including following columns:
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
#' # query gene of gencode version v26/hg38
#' geneInfo <- xQTLquery_gene("TP53")
#' geneInfo <- xQTLquery_gene(c("tp53","naDK","SDF4") )
#' geneInfo <- xQTLquery_gene(c("ENSG00000210195.2","ENSG00000078808"))
xQTLquery_gene <- function(genes="", geneType="auto", recordPerChunk=150){
  geneSymbol <- gencodeId <- entrezGeneId <- chromosome <- start <- end <- strand <- tss <- description <- cutF <- genesUpper <- NULL
  .<-NULL
  page_tmp <- 0
  pageSize_tmp <- recordPerChunk
  cutNum <- recordPerChunk
  genomeBuild="GRCh38/hg38"
  gencodeVersion="v26"

  # check genes
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }else if( any(duplicated(toupper(genes))) ){
    stop("Please remove duplicated genes.")
  }

  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("auto","geneSymbol", "gencodeId", "geneCategory")) ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\" or \"geneCategory\".")
  }

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else if( length(genes)==1 ){
      if( genes %in% gencodeGenetype[["V26"]] | genes %in% gencodeGenetype[["V19"]] ){
        geneType <- "geneCategory"
      }else{
        geneType <- "geneSymbol"
      }
    }else{
      geneType <- "geneSymbol"
    }
  }
  # message("geneType: ",geneType)

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

  message("Querying genes...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))

  ######################### if geneType is "geneCategory":
  if( geneType == "geneCategory" ){
    # Fetch all genes' info:
    gencodeGeneInfo <- xQTLquery_geneAll("v26")
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
      # bestFetchMethod <- apiAdmin_ping()
      # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
      #   return(NULL)
      # }
      # message("GTEx API successfully accessed!")
      # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])

      # download with "download" method and retry 3 times.
      url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

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
    }else if( length(genes)>1 ){
      # 分批下载：
      genesCut <- data.table::data.table(genes=genes, ID=1:length(genes), cutF = as.character(cut(1:length(genes),breaks=seq(0,length(genes)+cutNum,cutNum) )) )
      genesCut$genesUpper <- toupper(genesCut$genes)
      genesURL <- genesCut[,.(genesURL=paste0(genes,collapse = ",")),by=c("cutF")]
      if( any(unlist(lapply(genesURL$genesURL, nchar)) >3900) ){
        stop("Too many queried genes, please lower the value of \"recordPerChunk\", or reduce your input genes.")
      }
      outInfo <- data.table::data.table()

      # bestFetchMethod <- apiAdmin_ping()
      # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
      #   # message("Note: API server is busy or your network has latency, please try again later.")
      #   return(NULL)
      # }

      # message("GTEx API successfully accessed!")
      for(i in 1: nrow(genesURL) ){
        tmp_all <- data.table()
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

        # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
        # download with "download" method and retry 3 times.
        url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")


        url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
        if( nrow(url1GetText2Json2DT)==0 ){
          message( "0 record fatched!" )
          return(data.table::data.table())
        }
        url1GetText2Json2DT$genomeBuild <- genomeBuild
        tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
        tmp$genesUpper <- toupper(unlist(tmp[,geneType,with=FALSE]))
        tmp_all <- rbind(tmp_all, tmp)
        message("Batch ",i,"/",nrow(genesURL),". Downloaded  page", url1GetText2Json$page+1, "/",url1GetText2Json$numPages,"; ", length(na.omit(outInfo$gencodeId))+nrow(url1GetText2Json2DT), " records.")

        # if more pages:
        page_tmp<-page_tmp+1
        while( page_tmp <= (url1GetText2Json$numPages-1) ){
          url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                         "geneId=", genesURL[i,]$genesURL,"&",
                         "gencodeVersion=", gencodeVersion,"&",
                         "genomeBuild=",genomeBuild,"&",
                         "page=",page_tmp,"&",
                         "pageSize=", pageSize_tmp,"&",
                         "format=json"
          )
          url1 <- utils::URLencode(url1)
          # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
          url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
          url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
          if( nrow(url1GetText2Json2DT)==0 ){
            message( "0 record fatched!" )
            break()
          }
          url1GetText2Json2DT$genomeBuild <- genomeBuild
          tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
          tmp$genesUpper <- toupper(unlist(tmp[,geneType,with=FALSE]))
          tmp_all <- rbind(tmp_all, tmp)
          message("Batch ",i,". Downloaded  ", url1GetText2Json$page+1, "/",url1GetText2Json$numPages,"; ", length(na.omit(outInfo$gencodeId))+nrow(url1GetText2Json2DT), " records.")
          page_tmp <- page_tmp+1
        }
        # because of versioned and unversioned gencodeID, merge separately is needed!
        if(geneType == "gencodeId"){
          # versioned:
          tmp1 <- merge( genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp_all, by="genesUpper", sort = FALSE)[,-c("genesUpper")]
          # unversioned:
          tmp_all$genesUpper <- unlist(lapply(tmp_all$genesUpper, function(x){ stringr::str_split(x,stringr::fixed("."))[[1]][1] }))
          tmp2 <- merge( genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp_all, by="genesUpper", sort = FALSE)[,-c("genesUpper")]
          # combination:
          tmp_all <- merge(genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], rbind(tmp1, tmp2), by="genes",all.x=TRUE,sort = FALSE)
        }else{
          tmp_all <- merge( genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp_all, by="genesUpper",all.x=TRUE, sort = FALSE)[,-c("genesUpper")]
        }
        page_tmp=0
        outInfo <- rbind(outInfo, tmp_all)

        # message("Downloaded  ", round(i/nrow(genesURL)*100,2),"%; totally ", length(na.omit(outInfo$gencodeId)), " records fetched!")
        # rm(url1, url1Get, url1GetText, url1GetText2Json, url1GetText2Json2DT)
      }
      return(outInfo)
    }else{
      return(data.table::data.table())
    }
  }
}

#' @title Query sample's details with tissue name
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param dataType A character string. Options: "RNASEQ" (default), "WGS", "WES", "OMNI".
#' @param recordPerChunk (integer) number of records fetched per request (default: 200).
#' @param pathologyNotesCategories Default: pathologyNotes info is ignored.
#' @import data.table
#' @import utils
#' @import curl
#' @import jsonlite
#' @return returen a data.table object of samples' information
#' @export
#' @examples
#' sampleInfo <- xQTLquery_sampleByTissue("Brain - Amygdala" )
#' sampleInfo <- xQTLquery_sampleByTissue(tissueSiteDetail="Liver", pathologyNotesCategories=TRUE)
xQTLquery_sampleByTissue <- function( tissueSiteDetail="Liver", dataType="RNASEQ", recordPerChunk=200, pathologyNotesCategories=FALSE ){
  sampleId <- sex <- ageBracket <- pathologyNotes <- hardyScale <- NULL
  .<-NULL
  page_tmp <- 0
  pageSize_tmp <- as.integer(recordPerChunk)
  datasetId="gtex_v8"

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
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if( !(tissueSiteDetail %in% c("All", tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("0. All", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }
  # convert tissueSiteDetail to tissueSiteDetailId:
  tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId

  message("Querying samples by tissue...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))

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

  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])

  # download with "download" method and retry 3 times.
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

  if( ncol(url1GetText2Json$sample$pathologyNotesCategories)==0){
    pathologyNotesCategories <- FALSE
    message(" == No pathologyNotesCategories information found in tissue [", tissueSiteDetail, "] samples!" )
    tmp <- data.table::as.data.table( url1GetText2Json$sample[,which(names(url1GetText2Json$sample) != "pathologyNotesCategories")] )
    outInfo <- rbind(outInfo, tmp)
  }else{
    tmp <- data.table::as.data.table( url1GetText2Json$sample )
    outInfo <- rbind(outInfo, tmp)
  }

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

    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    # download with "download" method and retry 3 times.
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

    if( ncol(url1GetText2Json$sample$pathologyNotesCategories)==0){
      tmp <- data.table::as.data.table( url1GetText2Json$sample[,which(names(url1GetText2Json$sample) != "pathologyNotesCategories")] )
      outInfo <- rbind(outInfo, tmp)
    }else{
      tmp <- data.table::as.data.table( url1GetText2Json$sample )
      outInfo <- rbind(outInfo, tmp)
    }

    message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp <- page_tmp+1
  }
  if(pathologyNotesCategories){
    return(
      cbind(outInfo[,names(outInfo)[!str_detect(names(outInfo), stringr::regex("^pathologyNotesCategories."))], with=FALSE], outInfo[,names(outInfo)[str_detect(names(outInfo), stringr::regex("^pathologyNotesCategories."))], with=FALSE])
    )
  }else{
    return( outInfo[,names(outInfo)[!str_detect(names(outInfo), stringr::regex("^pathologyNotesCategories."))], with=FALSE] )
  }
}



#' @title Query sample's details with samples' IDs.
#' @param sampleIds A character vector or a string of sample ID.
#' @param recordPerChunk (integer) number of records fetched per request (default: 200).
#' @param pathologyNotesCategories Default: pathologyNotes info is ignored.
#'
#' @return a data.table object of samples' information.
#' @export
#'
#' @examples
#' sampleIds <- c("GTEX-11NUK-0011-R4a-SM-DO12B", "GTEX-11ONC-0011-R4b-SM-DO93H",
#'                "GTEX-11DXY-0526-SM-5EGGQ", "GTEX-13OVJ-1026-SM-5IFGI")
#' sampleInfo <- xQTLquery_sampleBySampleId(sampleIds)
xQTLquery_sampleBySampleId <- function(sampleIds,recordPerChunk=150, pathologyNotesCategories=FALSE ){
  . <- NULL

  sampleIds <- unique(sampleIds)
  if(!all(str_detect(sampleIds, "^GTEX-"))){
    stop("Samples ID should be start with \"GTEx-\", please check your input!")
  }

  pageSize_tmp <- as.integer(recordPerChunk)
  cutNum <- as.integer(recordPerChunk)

  # 分批下载防止样本量太大：
  samplesCut <- data.table::data.table(sampleIds=sampleIds, ID=1:length(sampleIds), cutF = as.character(cut(1:length(sampleIds),breaks=seq(0,length(sampleIds)+cutNum,cutNum) )) )
  samplesURL <- samplesCut[,.(samplesURL=paste0(sampleIds,collapse = ",")),by=c("cutF")]
  if( any(unlist(lapply(samplesURL$samplesURL, nchar)) >3900) ){
    stop("Too many queried genes, please lower the value of \"recordPerChunk\", or reduce your input samples")
  }


  message("Querying samples by sampleId...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))

  tmp_all <- data.table()
  for(i in 1: nrow(samplesURL) ){
    # 每个批次都得初始化page
    page_tmp <- 0
    message("== Batch ",i)

    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                   "sampleId=",samplesURL[i,]$samplesURL,"&",
                   "page=",page_tmp,"&",
                   "pageSize=", pageSize_tmp
    )
    url1 <- utils::URLencode(url1)
    outInfo <- data.table::data.table()

    # download with "download" method and retry 3 times.
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

    if( ncol(url1GetText2Json$sample$pathologyNotesCategories)==0){
      pathologyNotesCategories <- FALSE
      message(" == No pathologyNotesCategories information found." )
      tmp <- data.table::as.data.table( url1GetText2Json$sample[,which(names(url1GetText2Json$sample) != "pathologyNotesCategories")] )
      outInfo <- rbind(outInfo, tmp)
    }else{
      tmp <- data.table::as.data.table( url1GetText2Json$sample )
      outInfo <- rbind(outInfo, tmp)
    }
    message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp<-page_tmp+1

    while( page_tmp <= (url1GetText2Json$numPages-1)){
      url1 <- paste0("https://gtexportal.org/rest/v1/dataset/sample?",
                     "sampleId=",samplesURL[i,]$samplesURL,
                     "page=",page_tmp,"&",
                     "pageSize=", pageSize_tmp
      )
      url1 <- utils::URLencode(url1)

      # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
      # download with "download" method and retry 3 times.
      url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

      if( ncol(url1GetText2Json$sample$pathologyNotesCategories)==0){
        tmp <- data.table::as.data.table( url1GetText2Json$sample[,which(names(url1GetText2Json$sample) != "pathologyNotesCategories")] )
        outInfo <- rbind(outInfo, tmp)
      }else{
        tmp <- data.table::as.data.table( url1GetText2Json$sample )
        outInfo <- rbind(outInfo, tmp)
      }

      message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
      page_tmp <- page_tmp+1
    }
    if(pathologyNotesCategories){
      tmp_all <- rbind(tmp_all,cbind(outInfo[,names(outInfo)[!str_detect(names(outInfo), stringr::regex("^pathologyNotesCategories."))], with=FALSE], outInfo[,names(outInfo)[str_detect(names(outInfo), stringr::regex("^pathologyNotesCategories."))], with=FALSE]), fill=TRUE)
    }else{
      tmp_all<- rbind(tmp_all, outInfo[,names(outInfo)[!str_detect(names(outInfo), stringr::regex("^pathologyNotesCategories."))], with=FALSE])
    }
  }
  return(tmp_all)
}



#' @title Fetch details of all genes supported in GTEx.
#' @param gencodeVersion (character) options: "v26"(default, matched with gtex_v8) or "v19"
#' @param recordPerChunk (integer) number of records fetched per request (default: 2000).
#' @import utils
#' @import data.table
#' @return A data.table object of all genes' information.
xQTLquery_geneAll <- function(gencodeVersion="v26", recordPerChunk=2000){
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

  message("Querying all genes...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))

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
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # # message("GTEx API successfully accessed!")
  # suppressWarnings(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))

  # download with "download" method and retry 3 times.
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

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
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
    url1GetText2Json2DT$genomeBuild <- genomeBuild
    tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
    outInfo <- rbind(outInfo, tmp)
    message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp <- page_tmp+1
  }
  return(outInfo)
}


#' @title Query variant in GTEx with variant ID or dbSNP ID
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#'
#' @return A data.table object.
#' @export
#'
#' @examples
#' xQTLquery_varId("rs12596338")
#' xQTLquery_varId("chr11_66561248_T_C_b38")
xQTLquery_varId <- function(variantName="", variantType="auto"){
  datasetId="gtex_v8"
  ########## parameter check: variantName
  if(is.null(variantName) ||  any(is.na(variantName)) ){
    stop("Parameter \"variantName\" can not be NULL or NA!")
  }else if(length(variantName)!=1){
    stop("Length of \"variantName\" >1 !")
  }else if( variantType=="snpId" && !(stringr::str_detect(variantName, stringr::regex("^rs")))  ){
    stop("Parameter \"variantName\" must begin with a \"rs\", like: \"rs147502335\", \"rs147538909\".")
  }
  ########## parameter check: variantType
  if(is.null(variantType) ||  any(is.na(variantType)) ){
    stop("Parameter \"variantType\" can not be NULL or NA!")
  }else if(length(variantType)!=1){
    stop("Parameter \"variantType\" can not be NULL or NA!")
  }else if(  !(variantType %in% c("auto","snpId", "variantId"))  ){
    stop("Parameter \"variantType\" should be chosen from \"auto\", \"snpId\" and \"variantId\"!")
  }
  ########## parameter check: datasetId
  if(is.null(datasetId) ||  any(is.na(datasetId)) ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(length(datasetId)!=1){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }else if(  !(datasetId %in% c("gtex_v8", "gtex_v7"))  ){
    stop("Parameter \"datasetId\" should be chosen from \"gtex_v8\" or \"gtex_v7\"!")
  }

  # auto pick variantType
  if(variantType=="auto"){
    if(stringr::str_detect(variantName, stringr::regex("^rs"))){
      variantType <- "snpId"
    }else if(stringr::str_count(variantName,"_")==4){
      variantType <- "variantId"
    }else{
      stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
    }
  }

  ############# if variantType is "snpId":
  if( variantType == "snpId" ){
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/variant?format=json","&",
                   "datasetId=", datasetId,"&",
                   "snpId=", variantName)
    ############# if variantType is "variantId":
  }else if( variantType == "variantId" ){
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/variant?format=json","&",
                   "datasetId=", datasetId,"&",
                   "variantId=", variantName)
  }

  message("Querying variant...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))


  # url1 <- "https://gtexportal.org/rest/v1/dataset/variant?format=json&datasetId=gtex_v8&snpId=rs12596338"
  url1 <- utils::URLencode(url1)
  # test api server accessibility:
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # fetch data:
  # message("GTEx API successfully accessed!")
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2] )
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

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

#' @title Query varints in GTEx using genome position.
#' @param chrom (character) name of chromesome, including chr1-chr22, chrX, chrY.
#' @param pos An integer array.
#' @param recordPerChunk (integer) number of records fetched per request (default: 200).
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import utils
#'
#' @return A data.table object.
#' @export
#'
#' @examples
#' xQTLquery_varPos(chrom="chr1", pos=c(1102708,1105739))
xQTLquery_varPos <- function(chrom="", pos=numeric(0), recordPerChunk=200){
  .<-NULL
  datasetId="gtex_v8"
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

  message("Querying variant by position...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))


  var_tmp <- data.table::data.table(ID=1:length(pos), chrom=chrom, pos=pos)
  cutNum <- recordPerChunk
  ############### query with xQTLquery_varPos: START
  # var_tmp <- cbind(var_tmp, data.table::rbindlist(lapply(unique(outInfo$variantId), function(x){ splitOut = stringr::str_split(x,stringr::fixed("_"))[[1]];data.table::data.table(chrom=splitOut[1], pos=splitOut[2]) })))
  # query pos in batch per 100 terms.
  var_tmpCut <- cbind(var_tmp, data.table::data.table( ID=1:nrow(var_tmp), cutF = as.character(cut(1:nrow(var_tmp),breaks=seq(0,nrow(var_tmp)+cutNum,cutNum) )) ))
  var_tmpCut <- var_tmpCut[,.(posLis=list(pos)),by=c("chrom","cutF")]

  # check network:
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # out:
  outInfo <- data.table::data.table()
  for( ii in 1:nrow(var_tmpCut)){
    message("   Got varints: ",nrow(outInfo)+length(unlist(var_tmpCut[ii,]$posLis)),"/",nrow(var_tmp), " (",paste0(round( (nrow(outInfo)+length(unlist(var_tmpCut[ii,]$posLis)))/nrow(var_tmp),4)*100,"%"),")")
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/variant?format=json","&",
                   "datasetId=", datasetId,"&",
                   "chromosome=", chrom, "&",
                   "pos=",paste0(unlist(var_tmpCut[ii,]$posLis),collapse=",") )
    url1 <- utils::URLencode(url1)
    if(nchar(url1)>4000){
      stop("Too many positions were received, please reduce the number of positions or recordPerChunk.")
    }
    # fetch data:
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    tmp <- data.table::as.data.table(url1GetText2Json$variant)
    if(nrow(tmp)<1){
      message("No variants were found at thses positions.")
      return(NULL)
    }
    # eliminate the difference between the gtex_v7 and gtex_v8:
    if(datasetId == "gtex_v7"){
      tmp$b37VariantId <- tmp$variantId
    }
    tmp <- tmp[,c("variantId", "snpId","b37VariantId", "chromosome","pos", "ref", "alt","datasetId","maf01", "shorthand")]

    outInfo <- rbind(outInfo, tmp)
  }

  if(nrow(outInfo)==0){
    message("No variant found, please check your input.")
    return(data.table::data.table())
  }
  # print(paste("outInfo 行数：", nrow(outInfo)))

  outInfo <- outInfo[,c("variantId", "snpId","b37VariantId", "chromosome","pos", "ref", "alt","datasetId","maf01", "shorthand")]
  return(outInfo)
}


#' @title Fetch all details of a specified tissue or all tissues
#' @description
#'  Information includes tissue IDs, number of RNA-Seq samples, number of RNA-Seq samples with genotype, number of expressed genes, number of eGenes. Also includes tissueSiteDetail ID, name, abbreviation, uberon ID, and standard tissue colors. TissueSiteDetails are grouped by TissueSites. By default, this service reports from the latest GTEx release.
#' @param tissueName Tissue name, tissue ID or tissue site name. Default return all tissues' information. Can be choonse from `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#'
#' @return A data.table object.
#' @export
#'
#' @examples
#' tissueAll <- xQTLquery_tissue() # fetch all tissues
#' Brain <- xQTLquery_tissue("Brain")
#' Brains <- xQTLquery_tissue("Brain - Hippocampus")
xQTLquery_tissue <- function(tissueName=""){
  datasetId="gtex_v8"
  if(datasetId != "gtex_v8" && datasetId != "gtex_v7"){
    stop("\"datasetId\" must be choosen from \"gtex_v8\" or \"gtex_v7\"")
  }

  #  if tissue name is null:
  if(tissueName==""){
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/tissueInfo?format=json","&",
                   "datasetId=", datasetId)
    url1 <- utils::URLencode(url1)
    # bestFetchMethod <- apiAdmin_ping()
    # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    #   message("Note: API server is busy or your network has latency, please try again later.")
    #   return(NULL)
    # }
    # message("GTEx API successfully accessed!")
    # url1Get <- httr::GET(url1, httr::progress())
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

    url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$tissueInfo)
    return(url1GetText2Json2DT)
  }else{
    if(datasetId == "gtex_v8"){
      t1_tmp <- unique(na.omit(rbind(tissueSiteDetailGTExv8[tissueName,on="tissueSiteDetail"], tissueSiteDetailGTExv8[tissueName,on="tissueSiteDetailId"], tissueSiteDetailGTExv8[tissueName,on="tissueSite"])))
      if(nrow(t1_tmp)==0){
        stop("== [",tissueName, "] is not found in [",datasetId,"] please check your input.")
      }else{
        tissueSite = unique(t1_tmp$tissueSite)
      }
    }else if(datasetId == "gtex_v7"){
      t1_tmp <- unique(na.omit(rbind(tissueSiteDetailGTExv7[tissueName,on="tissueSiteDetail"], tissueSiteDetailGTExv7[tissueName,on="tissueSiteDetailId"],  tissueSiteDetailGTExv7[tissueName,on="tissueSite"])))
      if(nrow(t1_tmp)==0){
        stop("== [",tissueName, "] is not found in [",datasetId,"] please check your input.")
      }else{
        tissueSite = unique(t1_tmp$tissueSite)
      }
    }
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/tissueInfo?format=json","&",
                   "datasetId=", datasetId,"&",
                   "tissueSite=", tissueSite)
    url1 <- utils::URLencode(url1)
    # bestFetchMethod <- apiAdmin_ping()
    # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    #   message("Note: API server is busy or your network has latency, please try again later.")
    #   return(NULL)
    # }
    # message("GTEx API successfully accessed!")
    # url1Get <- httr::GET(url1, httr::progress())
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$tissueInfo)
    return(url1GetText2Json2DT)
  }
}


#' @title Heartbeat to check GTEx API server connectivity.
#' @description
#'  test GTEx API server and return download method.
#' @param fetchMethod fetchMethod.
#' @return A character string of fetchContent method.
#' @keywords internal
apiAdmin_ping <- function(fetchMethod=""){
  url1 <- "https://gtexportal.org/rest/v1/admin/ping"

  # platform-depend methods ordr:
  if( .Platform$OS.type =="unix"){
    methods_all <- c( "fromJSON", "download","curl")
  }else if(.Platform$OS.type =="windows"){
    methods_all <- c("download","curl",  "fromJSON")
  }else{
    methods_all <- c("download","curl",  "fromJSON")
  }


  if(fetchMethod==""){
    fetchMethod <- methods_all
  }else if( any(!fetchMethod %in% methods_all) ){
    stop("Invalid fetchMethod! ")
  }

  # download method:
  downloadMethod <-"auto"
  for( i in 1:length(fetchMethod)){
    tryCatch(
      {
        message("== Test GTEx network: ",fetchMethod[i])
        # if download:
        if( fetchMethod[i]=="download" ){

          if( !is.null(.Platform$OS.type) ){
            if( .Platform$OS.type =="windows"){ downloadMethod<-c("wininet", "libcurl" )}else{ downloadMethod<-  c("wget",  "libcurl"  )}
            if( .Platform$OS.type =="unix"){ downloadMethod<-c("wget", "libcurl")}
          }
          # start download:
            for(downM in 1:length(downloadMethod)){
              message( "    ==testing download method: ", downloadMethod[downM])
              if( downloadMethod[downM]=="libcurl" && !capabilities("libcurl") ){
                message(" Note: \"libcurl\" is not installed and can not be called ")
                next()
              }
              if( downloadMethod[downM]=="wget" && !capabilities("wget") ){
                stop(" Note: \"wget\" is not installed and can not be called, please install wget. ")
                next()
              }

              tryCatch({
                suppressWarnings(suppressMessages( outInfo <- fetchContent(url1, method = fetchMethod[i], downloadMethod = downloadMethod[downM]) ))
              },error=function(e){
                  # cat("ERROR :",conditionMessage(e), "\n")
                } )
              if( exists("outInfo") && outInfo=="Ping!"){
                # print(fetchMethod[i])
                message("== GTEx API successfully accessed!")
                return(c(fetchMethod[i], downloadMethod[downM]))
              }else{ next()}
            }
        }else{
          suppressWarnings(
            outInfo <- fetchContent(url1, method = fetchMethod[i])
          )
        }
        #
        if( exists("outInfo") && outInfo=="Ping!"){
          # print(fetchMethod[i])
          message("== GTEx API successfully accessed!")
          return(c(fetchMethod[i], downloadMethod))
        }else{
          next()
        }
      },
      # e = simpleError("test error"),
      error=function(cond){
        # message("   Method [",fetchMethod[i],"] failed!")
        return(NULL)
      },
      warning = function(cond){
        message(cond)
      }
    )
  }
  message("Note: GTEx API server is busy or your network has latency, please try again later.")
  return(NULL)
}

#' @title Heartbeat to check EBI API server connectivity.
#' @description
#'  test EBI API server and return download method.
#' @keywords internal
#' @return A character string of fetchContent method.
apiEbi_ping <- function(){
  url1 <- "https://www.ebi.ac.uk/eqtl/api/"
  fetchMethod = c("fromJSON","curl", "download")
  downloadMethod = "auto"
  for( i in 1:length(fetchMethod)){
    tryCatch(
      {
        message("== Test EBI network: ",fetchMethod[i])
        outInfo <- fetchContent(url1, method = fetchMethod[i], downloadMethod = downloadMethod)
        if( exists("outInfo") && !is.null(outInfo$`_links`) && length(outInfo$`_links`)>1 ){
          # print(fetchMethod[i])
          message("== EBI API successfully accessed!")
          return(c(fetchMethod[i], downloadMethod))
        }else{
          next()
        }
      },
      # e = simpleError("test error"),
      error=function(cond){
        # message("   Method [",fetchMethod[i],"] failed!")
        return(NULL)
      },
      warning = function(cond){
        message(cond)
      }
    )
  }
  message("Note: EBI API server is busy or your network has latency, please try again later.")
  return(NULL)
}


#' @title Fetch data using url by three methods
#'
#' @param url1 A url string.
#' @param method Can be chosen from "download", "curl", "fromJSON".
#' @param downloadMethod The same methods from utils::download.file function.
#' @param isJson Fetched content is a json file or not. Defaulst: TRUE.
#' @import utils
#' @import jsonlite
#' @import stringr
#' @import data.table
#' @importFrom curl curl_fetch_memory
#' @keywords internal
#' @return A json object.
fetchContent <- function(url1, method="curl", downloadMethod="auto", isJson=TRUE){
  # if( method == "GetWithHeader"){
  #   mycookie <- ''
  #   myheaders <- c('accept' ='ext/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8',
  #                  'accept-encoding' = 'gzip, deflate, br',
  #                  'accept-language' = 'zh-CN,zh;q=0.8,zh-TW;q=0.7,zh-HK;q=0.5,en-US;q=0.3,en;q=0.2',
  #                  'user-agent' = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:94.0) Gecko/20100101 Firefox/94.0',
  #                  'cookie' = mycookie)
  #   responseTmp <- httr::GET(url = url1, httr::add_headers(.headers = myheaders), encode="raw")
  #   #read_html()函数读入并解析网页内容
  #   webTmp <- rvest::read_html(responseTmp, encoding ="utf-8")
  #   url1GetText2Json<- jsonlite::fromJSON(rvest::html_text(webTmp))
  #   if( !exists("url1GetText2Json") || is.null(url1GetText2Json) || url1GetText2Json=="" || length(url1GetText2Json)==0 ){
  #     message("No data fetched, please check your input.")
  #     return(NULL)
  #   }
  #   return(url1GetText2Json)
  # }

  if(method == "fromJSON"){
    if(isJson){

      for (downloadTime in 1:4){
        # because of the limitation of the length of the url (<2084), so I rebulit the function.
        df <- try( url1GetText2Json <- jsonlite::fromJSON(url1, simplifyDataFrame=TRUE, flatten = TRUE), silent=TRUE)
        # methods::is(df, 'try-error')
        if( !(inherits(df, "try-error")) && exists("url1GetText2Json") ){
          break()
        }else if( (inherits(df, "try-error")) && !exists("url1GetText2Json") && downloadTime>3 ){
          message("No data fetched...")
          return(NULL)
        }
        if(downloadTime>1){message("=> Download failed and try again...",downloadTime-1,"/",3)}
        if(downloadTime>3){message("=> Your connection is unstable.")}
      }

      # url1GetText2Json <- jsonlite::fromJSON(url1, simplifyDataFrame=TRUE, flatten = TRUE)
      if( !exists("url1GetText2Json") || is.null(url1GetText2Json) || all(url1GetText2Json=="") || length(url1GetText2Json)==0 ){
        message("No data fetched, please check your input.")
        return(NULL)
      }
      return(url1GetText2Json)
    }else{
      message("fromJSON method return NULL!")
      return(NULL)
    }
  }else if( method == "rvest"){
    return(NULL)
    # if(isJson){
    #   url1GetText2Json <- jsonlite::fromJSON( rvest::html_text( rvest::read_html(url1) ) )
    #   if( !exists("url1GetText2Json") || is.null(url1GetText2Json) || all(url1GetText2Json=="") || length(url1GetText2Json)==0 ){
    #     message("No data fetched, please check your input.")
    #     return(NULL)
    #   }
    #   return(url1GetText2Json)
    # }else{
    #   url1GetText <- data.table::fread(rvest::html_text( rvest::read_html(url1) ),sep="\t", header = TRUE)
    # }

  }else if( method=="download" ){
    tmpFile <- tempfile(pattern = "file")
    if( file.exists(tmpFile) ){ file.remove(tmpFile) }
    # Retry for-loop R loop if error
    for (downloadTime in 1:4){
      if(downloadTime>1){message("=> Download failed and try again...",downloadTime-1,"/",3)}
      if(downloadTime>3){message("=> Your connection is unstable, please download  in brower directly using following url: ")
        message(url1)}
      df <- try(suppressWarnings(utils::download.file(url = url1, destfile=tmpFile, method=downloadMethod,quiet = TRUE )), silent=TRUE)
      if(!(inherits(df, "try-error"))) break
    }
    # suppressWarnings(utils::download.file(url = url1, destfile=tmpFile, method=downloadMethod,quiet = TRUE ))
    if(isJson){
      url1GetText2Json <-""
      if( file.exists(tmpFile) ){
        # url1GetText <- fread(tmpFile,sep="\n", header=FALSE)
        url1GetText <- readLines(tmpFile)
        # replace NAN with NULL:
        url1GetText <- stringr::str_replace_all(url1GetText, stringr::fixed("\":NaN,\""), "\":\"\",\"")
        url1GetText2Json <- jsonlite::fromJSON(url1GetText)
        if( !exists("url1GetText2Json") || is.null(url1GetText2Json) || all(url1GetText2Json=="") || length(url1GetText2Json)==0 ){
          message("No data fetched, please check your input.")
          return(NULL)
        }
        rm(url1GetText)
        return(url1GetText2Json)
      }else{
        stop("File download error.")
      }
    }else{
      url1GetText <- data.table::fread(tmpFile, sep="\t", header = TRUE)
      file.remove(tmpFile)
      # closeAllConnections()
      return(url1GetText)
    }
    close(file(tmpFile))
    if(file.exists(tmpFile)){ file.remove(tmpFile) }
  }else if( method == "curl"){
    url1Get <- curl::curl_fetch_memory(url1)
    url1GetText <- rawToChar(url1Get$content)
    # url1GetText <- RCurl::getURL(url1)
    # url1Get <- RCurl::getURLContent(url1)
    # if( url1Get$status_code!=200){
    #   mesage("Http status code: ", url1Get$status_code)
    # }

    if(isJson){
      # replace NAN with NULL:
      url1GetText <- stringr::str_replace_all(url1GetText, stringr::fixed("\":NaN,\""), "\":\"\",\"")
      url1GetText2Json <- jsonlite::fromJSON(url1GetText, flatten = FALSE)
      if( !exists("url1GetText2Json") || is.null(url1GetText2Json) || all(url1GetText2Json=="") || length(url1GetText2Json)==0 ){
        message("No data fetched, please check your input.")
        return(NULL)
      }
      return(url1GetText2Json)
    }else{
      url1GetText <- fread(url1GetText, sep="\t", header=TRUE)
      return(url1GetText)
    }
  }
}


#' @title Fetch records from
#'
#' @param url1 A url string.
#' @param method Can be chosen from "download", "curl".
#' @param downloadMethod The same methods from utils::download.file function.
#' @param termSize Number of records per request.
#' @param termStart Start position per request.
#' @import data.table
#' @keywords internal
#' @return A data.table object.
fetchContentEbi <- function(url1, method="fromJSON", downloadMethod="auto", termSize=1000, termStart=0){
  # method="curl"
  # downloadMethod="auto"
  # termSize=1000
  # termStart=0
  # i=1

  # message <- function(...){
  #   args <- list(...)
  #   argsTest <- paste0(args, collapse = "")
  #   stopifnot(require(crayon))
  #   warn <- magenta$underline
  #   cat( warn(argsTest),"\n" )
  # }

  embeddedList <- list()
  # If do not append size and start.
  if( termSize ==0 ){
    urlGot <- url1
    contentGot <- fetchContent(url1 = urlGot, method =method, downloadMethod = downloadMethod)
    embeddedList <- c(embeddedList,contentGot[[1]])
    embeddedListCount <- sum(unlist(lapply(embeddedList, function(x){ if(inherits(x, "data.frame")){return(nrow(x))}else{ return(length(x))} })))
    contentGotCount <- unlist(lapply(contentGot[[1]], function(x){ if(inherits(x, "data.frame")){return(nrow(x))}else{ return(length(x))} }))
    message("Got records: ",contentGotCount,"; Total records: ", embeddedListCount )

  # 这个API 有点问题，如果fetch 多个组织的结果，如果该组织剩下不足1000，就会将剩下的附加到这次fetch，导致 _links 里没有next，进而终止检索，所以再没有 next后，额外再fetch一次，如果有next则继续。
  }else if(termSize > 0){
    for(i in 1:99999999){
      urlGot <- paste0(url1, ifelse(stringr::str_detect(url1,stringr::fixed("?")),"&","?"),
                       "size=",as.integer(termSize),
                       "&start=",as.integer(termStart))
      contentGot <- fetchContent(url1 = urlGot, method =method, downloadMethod = downloadMethod)
      if( !is.null(contentGot$status) && contentGot$status==404 ){
        print(i)
        return(NULL)
      }else if( length(contentGot$`_embedded`$associations)==0 ){
        # is.null( contentGot$`_links` $`next`$href)
        embeddedList <- c(embeddedList,contentGot[[1]])
        embeddedListCount <- sum(unlist(lapply(embeddedList, function(x){ if(inherits(x, "data.frame")){return(nrow(x))}else{ return(length(x))} })))
        contentGotCount <- unlist(lapply(contentGot[[1]], function(x){ if(inherits(x, "data.frame")){return(nrow(x))}else{ return(length(x))} }))
        message("Request: ",i,"; Got records: ",contentGotCount,"; Total records: ", embeddedListCount)
        break()
      }else{
        embeddedList <- c(embeddedList,contentGot[[1]])
        embeddedListCount <- sum(unlist(lapply(embeddedList, function(x){ if(inherits(x, "data.frame")){return(nrow(x))}else{ return(length(x))} })))
        contentGotCount <- unlist(lapply(contentGot[[1]], function(x){ if(inherits(x, "data.frame")){return(nrow(x))}else{ return(length(x))} }))
        message("Request: ",i,"; Got records: ",contentGotCount,"; Total records: ", embeddedListCount)
        termStart <- termStart+as.integer(contentGotCount)
      }
    }
  }else{
    return(NULL)
  }
  if(length(embeddedList)>0){
    return(embeddedList)
  }else{
    return(NULL)
  }
}



#' @title retrieve snps from dbSNP using coordinate.
#'
#' @param chrom (character) name of chromesome, including chr1-chr22, chrX, chrY.
#' @param startPos A positive integer.
#' @param endPos A positive integer.
#' @param genomeBuild "GRCh38/hg38" or "GRCh38/hg19". Default: "GRCh38/hg38".
#' @param track "snp151Common", "snp150Common" or "snp147Common". Default: "snp151Common".
#' @import data.table
#' @import stringr
#' @import utils
#' @import curl
#' @import jsonlite
#' @keywords internal
#' @return A data.table object.
dbsnpQueryRange <- function(chrom="", startPos=-1, endPos=-1, genomeBuild="GRCh38/hg38", track="snp151Common" ){
  # url1<-"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=snp151Common;chrom=chr1;start=1;end=1000000"
  # chrom="chr1"
  # startPos=1
  # endPos=1000000
  # genomeBuild="GRCh38/hg38"
  # track="snp151Common"

  bestFetchMethod <- NULL

  # Parameter check: chrom
  if(  is.null(chrom) || is.na(chrom) || length(chrom)!=1 ){
    stop("Parameter \"chrom\" can not be NULL or NA!")
  }
  if( !stringr::str_detect(chrom, stringr::regex("^[cC][hH][rR]1[0-9]$|^[cC][hH][rR]2[0-3]$|^[cC][hH][rR][1-9]$|^[xXyY]$"))  ){
    stop("Parameter \"chrom\" stands for chromosome, which must be chosen from \"chr1-chr22, chrx, chry\" ")
  }

  # Parameter check: startPos
  if(  is.null(startPos) || is.na(startPos) || length(startPos)!=1 || !is.wholenumber(startPos) || startPos<0 ){
    stop("Parameter \"startPos\" must be a whole number.")
  }
  # Parameter check: endPos
  if(  is.null(endPos) || is.na(endPos) || length(endPos)!=1 || !is.wholenumber(endPos)|| endPos<0 ){
    stop("Parameter \"chrom\" must be a whole number.")
  }else{
    endPos <- as.integer(endPos)
  }
  # Parameter check: genomeBuild
  if(  is.null(genomeBuild) || is.na(genomeBuild) || length(genomeBuild)!=1 ){
    stop("Parameter \"chrom\" must be a whole number.")
  }else if(genomeBuild =="GRCh38/hg38"){
    genomeBuild="hg38"
  }else if(genomeBuild =="GRCh37/hg19"){
    genomeBuild="hg19"
  }else{
    stop("Parameter \"genomeBuild\" must be chosen form \"GRCh38/hg38\" and \"GRCh37/hg19\" ")
  }
  # Parameter check: track
  if(  is.null(track) || is.na(track) || length(track)!=1 ){
    stop("Parameter \"track\" can not be NULL or NA!")
  }else if( !track %in% c("snp151Common", "snp150Common", "snp147Common") ){
    stop("Parameter \"track\" should be chosen from \"snp151Common\", \"snp150Common\" and \"snp147Common\".")
  }

  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }

  # construct url:
  url1 <- paste0("https://api.genome.ucsc.edu/getData/track?",
                 "genome=",genomeBuild,
                 ";track=",track,
                 ";chrom=",chrom,
                 ";start=",startPos,
                 ";end=",endPos)

  url1 <- utils::URLencode(url1)
  message("  Downloading...")
  url1GetText2Json <- fetchContent(url1, method = "curl", downloadMethod = bestFetchMethod[2])
  outInfo <- data.table::as.data.table(url1GetText2Json[track][[track]])
  message("  Done")
  return(outInfo)
}


#' @title determine whether is a whole number:
#'
#' @param x A number
#' @param tol Don't change
#' @keywords internal
#' @return TRUE or FALSE
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}


#' @title EBIquery_allTerm
#'
#' @param term "associations", "molecular_phenotypes", "studies", "tissues", "qtl_groups", "genes" or "chromosomes".
#' @param termSize Number of fetched term.
#' @return A data.table object.
#' @export
#' @examples
#' \donttest{
#' # Fetch associatons:
#' associations <- data.table::rbindlist(EBIquery_allTerm("associations",termSize=0))
#'
#' # fetch molecular_phenotypes:
#' molecular_phenotypes <- EBIquery_allTerm("molecular_phenotypes", termSize=10)
#'
#' # fetch studies:
#' studies <- EBIquery_allTerm("studies")
#'
#' # fetch tissues:
#' tissues <- EBIquery_allTerm("tissues")
#'
#' # fetch tissue-study mapping relationships
#' tissue_S <- EBIquery_allTerm( paste0("tissues/", "UBER_0002046","/studies" ))
#'
#' # fetch qtl groups:
#' qtl_groups <- EBIquery_allTerm("qtl_groups")
#'
#' # Fetch genes:
#' geneList <- EBIquery_allTerm("genes", termSize=10)
#' }
EBIquery_allTerm <- function( term="genes", termSize=2000){
  bestFetchMethod <- apiEbi_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }else{
    message("EBI API server connected.")
  }

  message("Querying term...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))


  url1 <- "https://www.ebi.ac.uk/eqtl/api/"
  # allTerms <- fetchContent(url1,method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  # allTerms <- names(allTerms$`_links`)
  # if( !(term %in% allTerms) ){
  #   message("Parameter \"term\" must be chosen from \"", paste0(allTerms, collapse = "\", \""),".")
  # }
  #
  url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/", term)
  termInfo <- fetchContentEbi(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2], termSize = termSize)
  # studies:
  if(term == "studies"){
    termInfo <- data.table::rbindlist(lapply(termInfo$studies, function(x){ cbind(data.table(study_accession=x[[1]]),x[[2]]) }))
    termInfo <- termInfo[,c("study_accession")]
  }else if(term == "tissues"){
    termInfo <- termInfo$tissues
    termInfo$tissueType <- unlist(lapply(termInfo$tissue, function(x){ splitOut <- stringr::str_split(x, stringr::fixed("_"))[[1]]; splitOut[1] }))
  }else if(term =="genes"){
    termInfo <- do.call(rbind, lapply(termInfo, function(x){x['gene']}))
    row.names(termInfo) <- NULL
  }else{
    termInfo <- termInfo[[1]]
  }
  return(termInfo)
}



#' @title extract gene infor of specified genome from gencodeGeneInfoAllGranges
#'
#' @param gencodeGeneInfoAllGranges from internal data
#' @param genomeVersion "v26" (default) or "v19"
#' @import data.table
#' @import stringr
#' @importFrom BiocGenerics strand
#' @return A data.table object.
#' @export
#'
#' @examples
#' gencodeGeneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges)
extractGeneInfo <- function(gencodeGeneInfoAllGranges, genomeVersion="v26"){
  .<-NULL
  a <- data.table::copy(gencodeGeneInfoAllGranges)
  if(genomeVersion == "v26"){
    gencodeGeneInfo <- cbind(data.table::data.table(gencodeId=a$gencodeId, chromosome = as.character(seqnames(a)), strand=as.character(BiocGenerics::strand(a)) ), data.table::as.data.table(IRanges::ranges(a))[,.(start, end)])
    gencodeGeneInfo <- gencodeGeneInfo[start>0]
    # check:
    # nrow(gencodeGeneInfoV26[chromosome !="chrM"])
    # nrow(fintersect(gencodeGeneInfoV26[chromosome !="chrM",.(gencodeId=gencodeId_unversioned,chromosome, start,end, strand)], gencodeGeneInfo[,.(gencodeId, chromosome, start,end, strand)]))
  }else if( genomeVersion=="v19" ){
    gencodeGeneInfo <- cbind(data.table::data.table(gencodeId=a$gencodeId, chromosome = as.character(seqnames(a)), strand=as.character(BiocGenerics::strand(a)) ), data.table::as.data.table(IRanges::ranges(a$rangesV19))[,.(start, end)])
    gencodeGeneInfo <- gencodeGeneInfo[start>0]
    # check:
    # nrow(gencodeGeneInfoV19[chromosome !="chrMT"])
    # nrow(fintersect(gencodeGeneInfoV19[chromosome !="chrMT",.(gencodeId=gencodeId_unversioned,chromosome, start,end, strand)], gencodeGeneInfo[,.(gencodeId, chromosome, start,end, strand)]))
  }
  return(gencodeGeneInfo)
}



#' @title Retrive SNP pairwise LD from LDlink database
#'
#' @param targetSnp target SNP, support dbSNP IP.
#' @param population Supported population is consistent with the LDlink, which can be listed using function LDlinkR::list_pop()
#' @param windowSize Window around the highlighted snp for querying linkage disequilibrium information. Default:500000
#' @param method The same as fetchContent function, can be chosen from "download", "curl", "GetWithHeader".
#' @param genomeVersion "grch38"(default) or "grch37".
#' @param max_count To prevent download failure due to network fluctuations, max number of connection attempts.
#' @param token Ldlink token, default: "9246d2db7917"
#'
#' @export
#' @return A data.table object.
#'
#' @examples
#' # snpLD <- retrieveLD_LDproxy("rs3", windowSize=5000)
retrieveLD_LDproxy <- function(targetSnp="", population="EUR" , windowSize=50000, method="download", genomeVersion="grch38", max_count=3, token="9246d2db7917"){
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
    # message(url1)
    # url1 <- "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=rs3&pop=MXL&r2_d=r2&window=100000&genome_build=grch38&token=9246d2db7917"
    try( snpLDtmp <- fetchContent(url1, method=method, isJson=FALSE) )
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
