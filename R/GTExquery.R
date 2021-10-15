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
#'  # hg38:
#'  geneInfo <- apiRef_genes("TP53", "symbol", "v26", "GRCh38/hg38")
#'  geneInfo <- apiRef_genes(c("tp53","naDK","SDF4"), "symbol", "v26", "GRCh38/hg38")
#'  geneInfo <- apiRef_genes(c("ENSG00000210195.2","ENSG00000078808"),
#'                               geneType="gencodeId", "v26", "GRCh38/hg38")
#'  geneInfo <- apiRef_genes(c(51150,5590,4509), "entrezId", "v26", "GRCh38/hg38")
#'  # hg19:
#'  geneInfo <- apiRef_genes(c(51150,5590,4509), "entrezId", "v19","GRCh37/hg19")
#'  # get all protein_coding genes:
#'  proteinCoding <- apiRef_genes(genes="protein coding",
#'                                  geneType="geneCategory", "v26","GRCh38/hg38" )
#'  miRNA <- apiRef_genes(genes="miRNA", geneType="geneCategory", "v26","GRCh38/hg38" )
#'  }
apiRef_genes <- function(genes="", geneType="symbol", gencodeVersion="v26", genomeBuild="GRCh38/hg38"){
  # check null/na
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }
  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"symbol\", \"gencodeId\", \"entrezId\".")
  }else if( !(geneType %in% c("symbol", "gencodeId", "entrezId","geneCategory")) ){
    stop("Parameter \"geneType\" should be choosen from \"symbol\", \"gencodeId\", \"entrezId\",\"geneCategory\".")
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
    gencodeGeneInfo <- apiRef_geneAll("v26")
  }else if(gencodeVersion=="v19" & genomeBuild=="GRCh37/hg19"){
    # data(gencodeGeneInfoV19)
    # gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV19)
    gencodeGeneInfo <- apiRef_geneAll("v19")
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
    # geneCategory
  }else if(geneType == "geneCategory"){
    geneCategory = unique(gencodeGeneInfoV26$geneType)
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

#
#
#' @title Query gene info through API.
#'
#' @param genes gene symbol or gencode id (versioned or unversioned).
#' @param geneType
#' @param gencodeVersion "v26" or "v19"
#' @param recordPerChunk A integer value. Defaulut: 150
#'
#' @return A data.frame
#' @export
#'
#' @examples
#' GTExquery_gene(miRNA$geneSymbol,"v26", 150)
#'  miRNA <- apiRef_genes(genes="miRNA", geneType="geneCategory", "v26","GRCh38/hg38" )
#'  proT <- apiRef_genes(genes="protein coding", geneType="geneCategory", "v26","GRCh38/hg38" )
#'  proT1 <- GTExquery_gene(proT$gencodeId[1:1000],"symbol","v26", 150)
#'  proT2 <- GTExquery_gene(c("ENSG00000261459.1","ENSG00000283496.1"),"symbol","v26", 150)
GTExquery_gene <- function(genes="", geneType="symbol", gencodeVersion="v26", recordPerChunk=150){
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
    stop("Parameter \"geneType\" should be choosen from \"symbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("symbol", "gencodeId", "entrezId","geneCategory")) ){
    stop("Parameter \"geneType\" should be choosen from \"symbol\", \"gencodeId\".")
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

  # API处理GENE，单个字符串模糊匹配，多个字符串精确匹配:
  if(length(genes)==1){
    # construct url:
    url1 <- paste0("https://gtexportal.org/rest/v1/reference/gene?",
                   "geneId=^", genes,"$&",
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
      outInfo <- cbind(data.table(genes=genes), outInfo)
      message( nrow(tmp), " record has been obtained!" )
      return(outInfo)
    }else{
      message("")
      return(data.table::data.table())
    }
  }else if( length(genes)>1 ){
    # 分批下载：
    genesCut <- data.table(genes=genes, ID=1:length(genes), cutF = as.character(cut(1:length(genes),breaks=seq(0,length(genes)+cutNum,cutNum) )) )
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
      url1GetText2Json2DT$genomeBuild <- genomeBuild
      tmp <- url1GetText2Json2DT[,.(geneSymbol, gencodeId, entrezGeneId, geneType, chromosome, start, end, strand, tss, gencodeVersion,genomeBuild, description)]
      tmp$genesUpper <- toupper(unlist(tmp[,geneType,with=FALSE]))
      # rm(url1, url1Get, url1GetText, url1GetText2Json, url1GetText2Json2DT)
      tmp <- merge(genesCut[cutF==genesURL[i,]$cutF,.(genes, genesUpper)], tmp, by="genesUpper",all.x=TRUE, sort = FALSE)[,-c("genesUpper")]
      outInfo <- rbind(outInfo, tmp)
      message("Downloaded part ", i, "/",nrow(genesURL),"; ", nrow(outInfo), " records.")
    }
    return(outInfo)
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
#'   #
#'   GTExquery_sample( tissueSiteDetail="Liver", dataType="RNASEQ", datasetId="gtex_v8", pageSize=200 )
#'   GTExquery_sample( tissueSiteDetail="All", dataType="RNASEQ", datasetId="gtex_v8",pageSize=2000 )
#'   GTExquery_sample( "Adipose - Visceral (Omentum)", "RNASEQ", "gtex_v8", 200 )
#'   }
GTExquery_sample <- function( tissueSiteDetail="Liver", dataType="RNASEQ", datasetId="gtex_v8", pageSize=200 ){
  page_tmp <- 0
  pageSize_tmp <- as.integer(pageSize)

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
#' @param toSummarizedExperiment return a data.frame or a summarizedExperiment object. Default return a data.frame.
#' @import data.table
#' @import curl
#' @import stringr
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom   IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @return return a SummarizedExperiment or data.frame object harbing gene expression profiles.
#' @export
#' @examples
#' \donttest{
#'   expProfiles <- GTExquery_exp(c("ENSG00000210195.2","ENSG00000078808"),
#'                                  "gencodeId", "Liver", "gtex_v8")
#'   expProfiles <- GTExquery_exp(c("tp53","naDK","SDF4"),
#'                                  "symbol", "Artery - Coronary", "gtex_v8",
#'                                  toSummarizedExperiment=TRUE)
#'   expProfiles <- GTExquery_exp(c("tp53","naDK","SDF4"),
#'                                  "symbol", "Artery - Coronary", "gtex_v7")
#'   # get miRNA expression:
#'   miRNA <- apiRef_genes (genes="miRNA", geneType="geneCategory", "v26","GRCh38/hg38" )
#'   miRNAExp <- GTExquery_exp(miRNA$gencodeId[1:10], geneType = "gencodeId")
#'   proT <- apiRef_genes (genes="protein coding", geneType="geneCategory", "v26","GRCh38/hg38" )
#'   proTexp <- GTExquery_exp(proT$gencodeId[1:10], geneType = "gencodeId")
#'   }
GTExquery_exp <- function(genes="", geneType="symbol", tissueSiteDetail="Liver", datasetId="gtex_v8", toSummarizedExperiment=TRUE, recordPerChunk=150  ){
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

  # convert genes. parameter check is unnecessary for this, because GTExquery_gene check it internally.
  message("== Querying gene info from API server:")
  geneInfo <- GTExquery_gene(genes=genes, geneType=geneType, gencodeVersion=gencodeVersion, recordPerChunk=recordPerChunk)
  message("== Done.")
  # geneInfo <-  GTExquery_gene(c("tp53","naDK","SDF4"), "symbol", "v26", "GRCh38/hg38")

  # get sample info:
  message("== Fetching sample information from API server:")
  sampleInfo <- GTExquery_sample(tissueSiteDetail=tissueSiteDetail, dataType="RNASEQ", datasetId=datasetId, pageSize=recordPerChunk )
  message("== Done.")
  if( is.null(sampleInfo)||!exists("sampleInfo") ){
    stop("Failed to fetch sample information.")
  }

  ############### 分批下载：
  geneInfo <- cbind(geneInfo, data.table(cutF = cut(1:nrow(geneInfo),breaks=seq(0,nrow(geneInfo)+cutNum,cutNum) ) ))
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


  ################# construct url:
  url1 <- paste0("https://gtexportal.org/rest/v1/expression/geneExpression?",
                 "datasetId=", datasetId,"&",
                 "gencodeId=", paste(geneInfo$gencodeId, collapse = ","), "&",
                 "tissueSiteDetailId=", tissueSiteDetailId,"&",
                 "&format=json"
  )
  url1 <- utils::URLencode(url1)
  message("== Fetching sample information from API server:")
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
    if(nrow(outInfo)==0){
      message("No expression profiles were found in ",tissueSiteDetail, " of thess ", length(genes), " genes!")
      message("== Done.")
      return(data.table::data.table())
    }
    outInfo <- merge(geneInfo, outInfo, by= c("geneSymbol", "gencodeId"))
    outInfo <- cbind(outInfo[,c("genes")], outInfo[,-c("genes")])
    expDT <- data.table::as.data.table(t(data.table::as.data.table(outInfo$data)))
    data.table::setDF(expDT)
    # some genes were excluded from gtex:
    excludedGenes <- data.table::fsetdiff(geneInfo[,.(gencodeId)],outInfo[,.(gencodeId)])
    includedGenes <- data.table::fintersect(geneInfo[,.(gencodeId)],outInfo[,.(gencodeId)])
    includedGenes <- merge(geneInfo, includedGenes, by="gencodeId", sort=FALSE)
    # note: includedGenes replace geneInfo in following:
    rownames(expDT) <- includedGenes$geneSymbol
    colnames(expDT) <- sampleInfo$sampleId
    message( ncol(expDT)," Expression profiles of ",length(unique(outInfo$gencodeId)), " genes in ",unique(outInfo$tissueSiteDetailId)," were obtained."," Unit: ",paste0(unique(outInfo$unit), collapse = ","),"." )
    #
    if( toSummarizedExperiment ){      # construc summarizedExperiment object::
      # gene info:
      expRowRanges <- GenomicRanges::GRanges(includedGenes$chromosome,
                                             IRanges::IRanges(includedGenes$start, includedGenes$end),
                                             strand= includedGenes$strand,
                                             includedGenes[,.(genes, geneSymbol, gencodeId, entrezGeneId, geneType, tss, description)]
                                             # ,seqinfo = GenomeInfoDb::Seqinfo(genome= stringr::str_split(unique(includedGenes$genomeBuild)[1],stringr::fixed("/"))[[1]][2])
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
      message("== Done.")
      return(outInfoSE)
    }else{
      outInfoDF <- cbind(outInfo[,-c("data")], expDT)
      data.table::setDF(outInfoDF)
      message("== Done.")
      return(outInfoDF)
    }
  }else{
    message("")
    return(data.frame())
  }
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

# Indepedently fetch geneInfo:
#' Title
#'
#' @param gencodeVersion A character string. "v26" or "v19".
#' @param recordPerChunk A integer value. From 1 to 2000.
#'
#' @return
#' @import utils
#' @import data.table
#' @export
#'
#' @examples
#' /donttest{
#'  gencodeGeneInfoV26 <- apiRef_geneAll("v26")
#'  gencodeGeneInfoV19 <- apiRef_geneAll("v19")
#' }
apiRef_geneAll <- function(gencodeVersion="v26", recordPerChunk=2000){
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
    message("")
    return(data.table::data.table())
  }
}









