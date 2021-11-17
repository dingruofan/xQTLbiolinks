#' #' @title coloc with a gwas dataset
#' #'
#' #' @param gwasDF
#' #' @param traitGene
#' #' @param traitGeneType
#' #' @param tissueSiteDetail
#' #' @param study_id
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' #' \donttest{
#' #'  gwasDF_raw <- data.table::fread("D:\\R_project\\GWAS_Type-2-Diabetes_Wood_2016.txt.gz", sep="\t", header=TRUE)
#' #'  # gwasDF <- gwasDF_raw[,.(rsid, pvalue, beta, varbeta)]
#' #'  names(gwasDF) <- c("snpId", "pValue", "beta", "varbeta")
#' #'  gwasDF <- gwasDF_raw[,.(rsid, pvalue, Freq)]
#' #'  names(gwasDF) <- c("snpId", "pValue", "freq")
#' #' }
#' GTExanalyze_coloc <- function(gwasDF, traitGene="", traitGeneType="geneSymbol", tissueSiteDetail="", study_id="gtex_v8"){
#'   # traitGene = "CAMK1D"
#'   # traitGeneType="geneSymbol"
#'   # tissueSiteDetail="Lung"
#'   # study_id="gtex_v8"
#'
#'   pos <- pValue.gwas <- pValue.etql <-NULL
#'   gencodeVersion <- "v26"
#'   if( all(c("snpId", "pValue", "freq") %in% names(gwasDF)) || all(c("snpId", "pValue", "beta", "varbeta") %in% names(gwasDF)) ){
#'
#'   }else{
#'
#'   }
#'
#'   # parameter check: GWAS snpId
#'   message("==Checking parameter!")
#'   # random check:
#'   if( nrow(gwasDF)>10000){
#'     set.seed(1)
#'     randomId <- sample(nrow(gwasDF),10000,replace = FALSE)
#'     if( !all(unlist(lapply(gwasDF$snpId[randomId], function(x){ stringr::str_detect(x,stringr::regex("^rs[0-9]{1,30}[0-9]$")) }))) ){
#'       message("The first column of \"gwasDF\" must be snp ID, which start with \"rs\"! ")
#'       return(data.table::data.table())
#'     }
#'     rm(randomId)
#'   }
#'
#'   if( is.null(traitGeneType) ||  any(is.na(traitGeneType)) || any(traitGeneType=="") || length(traitGeneType)!=1 ){
#'     stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
#'   }else if( !(traitGeneType %in% c("geneSymbol", "gencodeId")) ){
#'     stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
#'   }
#'   message("== Done")
#'
#'   #
#'   message("== Querying gene info from API server:")
#'   # check network:
#'   bestFetchMethod <- apiAdmin_ping()
#'   if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
#'     # message("Note: API server is busy or your network has latency, please try again later.")
#'     return(NULL)
#'   }
#'   suppressMessages( geneInfo <- GTExquery_gene(traitGene, geneType = traitGeneType, gencodeVersion = gencodeVersion) )
#'   if( exists("geneInfo") && nrow(geneInfo)>0 ){
#'     message("== Done")
#'   }else if( nrow(geneInfo)>1 ){
#'     message("Overall ",nrow(geneInfo)," genes were detected: [",paste0(geneInfo$gencodeId, collapse = ","),"]. ", "Please select one of the genes and rerun the function.")
#'   }else{
#'     message("Can not get gene info for [", traitGene,"] in ", study_id,".")
#'   }
#'
#'   ######## Fetch all associations:
#'   eqtlAsso <- GTExdownload_assoAll(gene=geneInfo$gencodeId, geneType = "gencodeId", tissueSiteDetail = tissueSiteDetail )
#'   # if no eqtl got:
#'   if( nrow(eqtlAsso)==0 || !exists("eqtlAsso")){
#'     stop("No eqtl associations were found for ",traitGeneType,": [",traitGene,"].")
#'   }
#'
#'   # merge and obtain info：
#'   gwas_eqtl <- merge(gwasDF, eqtlAsso, by="snpId", sort=FALSE, suffixes = c(".gwas",".etql"))
#'   if( nrow(gwas_eqtl)==0 || !exists("gwas_eqtl") ){
#'     message("No intersection of GWAS and eQTL dataset!")
#'     return(data.table::data.table())
#'   }
#'
#'   result <- coloc.abf(dataset1=list(pvalues=gwas_eqtl$pValue.gwas, type="cc", s=0.33, N=50000),
#'                       dataset2=list(pvalues=gwas_eqtl$pValue.etql, type="quant", N=10000), MAF=gwas_eqtl$maf)
#'
#'
#' }
