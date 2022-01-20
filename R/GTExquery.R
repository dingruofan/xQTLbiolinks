#' @title Query gene information through API.
#' @description
#'  users can use this function to search gene of interest
#' @param genes A charater vector or a string of gene symbol, gencode id (versioned), or a charater string of gene type.
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
#' @param geneType A character string. "geneSymbol"(default), "gencodeId" or "geneCategory".
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
  }else if( any(duplicated(genes)) ){
    stop("Please remove duplicated genes.")
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
      bestFetchMethod <- apiAdmin_ping()
      if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
        return(NULL)
      }
      # message("GTEx API successfully accessed!")
      url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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
      bestFetchMethod <- apiAdmin_ping()
      if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
        # message("Note: API server is busy or your network has latency, please try again later.")
        return(NULL)
      }
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
        url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
        url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
        if( nrow(url1GetText2Json2DT)==0 ){
          message( "0 record fatched!" )
          return(data.table::data.table())
        }
        url1GetText2Json2DT$genomeBuild <- genomeBuild
        tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
        tmp$genesUpper <- toupper(unlist(tmp[,geneType,with=FALSE]))
        tmp_all <- rbind(tmp_all, tmp)
        message("Batch ",i,". Downloaded  ", url1GetText2Json$page+1, "/",url1GetText2Json$numPages,"; ", length(na.omit(outInfo$gencodeId))+nrow(url1GetText2Json2DT), " records.")

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
          url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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
  }else if(length(tissueSiteDetail)!=1){
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
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  # message("GTEx API successfully accessed!")
  url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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
    url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    tmp <- data.table::as.data.table(url1GetText2Json$sample)
    outInfo <- rbind(outInfo, tmp)
    message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp <- page_tmp+1
  }
  return(outInfo[,.(sampleId, sex, ageBracket,datasetId, tissueSiteDetail, tissueSiteDetailId, pathologyNotes, hardyScale,  dataType )])
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
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  # message("GTEx API successfully accessed!")
  suppressWarnings(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
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
    url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$gene)
    url1GetText2Json2DT$genomeBuild <- genomeBuild
    tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
    outInfo <- rbind(outInfo, tmp)
    message("Total records: ",url1GetText2Json$recordsFiltered,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$numPages)
    page_tmp <- page_tmp+1
  }
  return(outInfo)
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
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  # fetch data:
  # message("GTEx API successfully accessed!")
  url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2] )
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
#' @param recordPerChunk A integer value (1-200). number of records fetched per request (default: 200).
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import utils
#'
#' @return A data.table.
#' @export
#'
#' @examples
#' \donttest{
#'  GTExquery_varPos(chrom="chr1", pos=c(1102708,1105739),"gtex_v8")
#'  GTExquery_varPos(chrom="1", pos=c(1038088,1041119),"gtex_v7")
#'  GTExquery_varPos("1", c(1246438, 1211944, 1148100),"gtex_v7")
#' }
GTExquery_varPos <- function(chrom="", pos=numeric(0), datasetId="gtex_v8", recordPerChunk=200){
  .<-NULL
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

  var_tmp <- data.table::data.table(ID=1:length(pos), chrom=chrom, pos=pos)
  cutNum <- recordPerChunk
  ############### query with GTExquery_varPos: START
  # var_tmp <- cbind(var_tmp, data.table::rbindlist(lapply(unique(outInfo$variantId), function(x){ splitOut = stringr::str_split(x,stringr::fixed("_"))[[1]];data.table::data.table(chrom=splitOut[1], pos=splitOut[2]) })))
  # query pos in batch per 100 terms.
  var_tmpCut <- cbind(var_tmp, data.table::data.table( ID=1:nrow(var_tmp), cutF = as.character(cut(1:nrow(var_tmp),breaks=seq(0,nrow(var_tmp)+cutNum,cutNum) )) ))
  var_tmpCut <- var_tmpCut[,.(posLis=list(pos)),by=c("chrom","cutF")]

  # check network:
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  # out:
  outInfo <- data.table::data.table()
  for( ii in 1:nrow(var_tmpCut)){
    message("   Got varints: ",nrow(outInfo)+length(unlist(var_tmpCut[ii,]$posLis)),"/",nrow(var_tmp))
    url1 <- paste0("https://gtexportal.org/rest/v1/dataset/variant?format=json","&",
                   "datasetId=", datasetId,"&",
                   "chromosome=", chrom, "&",
                   "pos=",paste0(unlist(var_tmpCut[ii,]$posLis),collapse=",") )
    url1 <- utils::URLencode(url1)
    if(nchar(url1)>4000){
      stop("Too many positions were received, please reduce the number of positions or recordPerChunk.")
    }
    # fetch data:
    url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
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



#' @title Heartbeat to check GTEx API server connectivity.
#' @description
#'  test GTEx API server and return download method.
#' @param fetchMethod fetchMethod.
#' @export
#' @return A character of fetchContent method.
#' @examples
#' \donttest{
#'   apiAdmin_ping()
#'  }
apiAdmin_ping <- function(fetchMethod=""){
  url1 <- "https://gtexportal.org/rest/v1/admin/ping"
  methods_all <- c("download","curl", "GET", "fromJSON")
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
          if( !capabilities("libcurl")){
            message(" \"libcurl\" can not be called, please install curl package  ")
            return(NULL)}
          if( !is.null(.Platform$OS.type) ){
            if( .Platform$OS.type =="windows"){ downloadMethod<-c("wininet", "libcurl" )}else{ downloadMethod<-  c("wget",  "libcurl"  )}
          }
          # start download:
            for(downM in 1:length(downloadMethod)){
              message( "    ==testing download method: ", downloadMethod[downM])
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
#' @export
#' @return A character of fetchContent method.
#' @examples
#' \donttest{
#'   apiEbi_ping()
#'  }
apiEbi_ping <- function(){
  url1 <- "https://www.ebi.ac.uk/eqtl/api/"
  fetchMethod = c("fromJSON","curl", "download","GET")
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
#' @param method Can be chosen from "download", "curl", "GetWithHeader", or "GET".
#' @param downloadMethod The same methods from utils::download.file function.
#' @param isJson Fetched content is a json file or not. Defaulst: TRUE.
#' @import utils
#' @import jsonlite
#' @import stringr
#' @import data.table
#' @importFrom curl curl_fetch_memory
#' @importFrom httr GET
#' @return A json object.
#' @export
#'
#' @examples
#' \donttest{
#'  url1 <- "https://gtexportal.org/rest/v1/admin/ping"
#'  fetchContent(url1, method="download")
#'  url1 <- paste0("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?",
#'                 "var=rs3&pop=MXL&r2_d=r2&window=500000&genome_build=grch38&token=9246d2db7917")
#'  fetchContent(url1, method="download", isJson=FALSE)
#' }
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
      url1GetText2Json <- jsonlite::fromJSON(url1, simplifyDataFrame=TRUE, flatten = TRUE)
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
    utils::download.file(url = url1, destfile=tmpFile, method=downloadMethod,quiet = TRUE )
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
        file.remove(tmpFile)
        return(url1GetText2Json)
      }else{
        stop("File download error.")
      }
    }else{
      url1GetText <- data.table::fread(tmpFile, sep="\t", header = TRUE)
      return(url1GetText)
    }
  }else if( method == "curl"){
    url1Get <- curl::curl_fetch_memory(url1)
    # if( url1Get$status_code!=200){
    #   mesage("Http status code: ", url1Get$status_code)
    # }
    url1GetText <- rawToChar(url1Get$content)
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
  }else if( method == "GET" ){
    url1Get <- httr::GET(url1)
    # if( url1Get$status_code!=200){
    #   message("Http status code: ", url1Get$status_code)
    # }
    url1GetText <- rawToChar(url1Get$content)
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
  }else{
    stop("please choose the right download method.")
  }
}


#' @title Fetch records from
#'
#' @param url1 A url string.
#' @param method Can be chosen from "download", "curl" or "GET".
#' @param downloadMethod The same methods from utils::download.file function.
#' @param termSize Number of records per request.
#' @param termStart Start position per request.
#' @import data.table
#' @return A data.frame
#' @export
#'
#' @examples
#' \donttest{
#'  url1 <- "https://www.ebi.ac.uk/eqtl/api/tissues/CL_0000057/associations?gene_id=ENSG00000141510"
#'  gtexAsoo <- fetchContentEbi(url1)
#' }
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
  # Don't append size and start.
  if( termSize ==0 ){
    urlGot <- url1
    contentGot <- fetchContent(url1 = urlGot, method =method, downloadMethod = downloadMethod)
    embeddedList <- c(embeddedList,contentGot[[1]])
    embeddedListCount <- sum(unlist(lapply(embeddedList, function(x){ if(class(x)=="data.frame"){return(nrow(x))}else{ return(length(x))} })))
    contentGotCount <- unlist(lapply(contentGot[[1]], function(x){ if(class(x)=="data.frame"){return(nrow(x))}else{ return(length(x))} }))
    message("Got records: ",contentGotCount,"; Total records: ", embeddedListCount )
  }else if(termSize > 0){
    for(i in 1:1000000){
      urlGot <- paste0(url1, ifelse(stringr::str_detect(url1,stringr::fixed("?")),"&","?"),
                       "size=",as.integer(termSize),
                       "&start=",as.integer(termStart))
      contentGot <- fetchContent(url1 = urlGot, method =method, downloadMethod = downloadMethod)
      if( !is.null(contentGot$status) && contentGot$status==404 ){
        print(i)
        return(NULL)
      }else if( is.null(contentGot$`_links`$`next`$href) ){
      # }else if( length(contentGot[[1]][[1]])==0|| length(contentGot[[1]][[1]][[1]])==0 ){
        embeddedList <- c(embeddedList,contentGot[[1]])
        embeddedListCount <- sum(unlist(lapply(embeddedList, function(x){ if(class(x)=="data.frame"){return(nrow(x))}else{ return(length(x))} })))
        contentGotCount <- unlist(lapply(contentGot[[1]], function(x){ if(class(x)=="data.frame"){return(nrow(x))}else{ return(length(x))} }))
        message("Requset: ",i,"; Got records: ",contentGotCount,"; Total records: ", embeddedListCount)
        break()
      }else{
        embeddedList <- c(embeddedList,contentGot[[1]])
        embeddedListCount <- sum(unlist(lapply(embeddedList, function(x){ if(class(x)=="data.frame"){return(nrow(x))}else{ return(length(x))} })))
        contentGotCount <- unlist(lapply(contentGot[[1]], function(x){ if(class(x)=="data.frame"){return(nrow(x))}else{ return(length(x))} }))
        message("Requset: ",i,"; Got records: ",contentGotCount,"; Total records: ", embeddedListCount)
        termStart <- termStart+as.integer(termSize)
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
#' @param chrom A character string. Chromesome of human, including chr1-chr23, chrX, chrY.
#' @param startPos A positive integer.
#' @param endPos A positive integer.
#' @param genomeBuild "GRCh38/hg38" or "GRCh38/hg19". Default: "GRCh38/hg38".
#' @param track "snp151Common", "snp150Common" or "snp147Common". Default: "snp151Common".
#' @import data.table
#' @import stringr
#' @import utils
#' @import curl
#' @import jsonlite
#' @return A data.table.
#' @export
#'
#' @examples
#' \donttest{
#'  snpInfo <- dbsnpQueryRange(chrom="chr1", startPos=1,
#'    endPos=10000000, genomeBuild="GRCh38/hg38", track="snp151Common" )
#' }
dbsnpQueryRange <- function(chrom="", startPos=-1, endPos=-1, genomeBuild="GRCh38/hg38", track="snp151Common" ){
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
    stop("Parameter \"chrom\" stands for chromosome, which must be chosen from \"chr1-chr23, chrx, chry\" ")
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


#' @title determine is a whole number:
#'
#' @param x A number
#' @param tol Don't change
#' @export
#' @return TRUE or FALSE
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}


#' @title EBIquery_allTerm
#'
#' @param term "associations", "molecular_phenotypes", "studies", "tissues", "qtl_groups", "genes" or "chromosomes".
#' @param termSize Number of fetched term.
#' @return a data.table
#' @export
#'
#' @examples
#' \donttest{
#'
#'  associations <- data.table::rbindlist(EBIquery_allTerm("associations",termSize=0))
#'  # molecular_phenotypes <- EBIquery_allTerm("molecular_phenotypes")
#'  studies <- EBIquery_allTerm("studies")
#'  tissues <- EBIquery_allTerm("tissues")
#'  qtl_groups <- EBIquery_allTerm("qtl_groups")
#'  # genes <- EBIquery_allTerm("genes")
#'  # chromosomes <- EBIquery_allTerm("chromosomes")
#'
#'  merge(qtl_groups, tissueSiteDetailGTExv8, by.x="qtl_group", by.y="tissueSiteDetail")
#' }
EBIquery_allTerm <- function( term="genes",termSize=5000){
  bestFetchMethod <- apiEbi_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }else{
    message("EBI API server connected.")
  }
  url1 <- "https://www.ebi.ac.uk/eqtl/api/"
  allTerms <- fetchContent(url1,method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  allTerms <- names(allTerms$`_links`)
  if( !(term %in% allTerms) ){
    message("Parameter \"term\" must be chosen from \"", paste0(allTerms, collapse = "\", \""),".")
  }
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
#' @return a data.table
#' @export
#'
#' @examples
#' \donttest{
#'   gencodeGeneInfo <- extractGeneInfo(gencodeGeneInfoAllGranges)
#' }
extractGeneInfo <- function(gencodeGeneInfoAllGranges, genomeVersion="v26"){
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




