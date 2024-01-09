#' @title Download normalized gene expression at the sample level for a specified tissue.
#' @param genes (character string or a character vector) gene symbols or gencode ids (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param toSummarizedExperiment a logical value indicating whether to return a data.frame or a summarizedExperiment object. Default: TRUE, return a toSummarizedExperiment object.
#' @param recordPerChunk (integer) number of records fetched per request (default: 80).
#' @param pathologyNotesCategories a logical value indicating whether to return pathologyNotes. Default: FALSE, the pathologyNotes is ignored.
#' @import data.table
#' @import stringr
#' @import jsonlite
#' @import stats
#' @importFrom GenomicRanges GRanges
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom   IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @return return a SummarizedExperiment or a data.table object harbing gene expression profiles and samples' information.
#' @export
#' @examples
#' \donttest{
#' # Download gene expression with a genecode ID:
#' expProfiles <- xQTLdownload_exp("ENSG00000210195.2", tissueSiteDetail="Liver")
#'
#' # Download gene expression into a SummarizedExperiment object:
#' expProfiles <- xQTLdownload_exp("ENSG00000210195.2", tissueSiteDetail="Liver",
#'                toSummarizedExperiment=TRUE)
#' # extract expression profile from SummarizedExperiment object:
#' expDT <- SummarizedExperiment::assay(expProfiles)
#' # extract samples' detail from SummarizedExperiment object:
#' sampleDT <- SummarizedExperiment::colData(expProfiles)
#'
#' # Download gene expression profiles of multiple genes:
#' expProfiles <- xQTLdownload_exp(c("tp53","naDK","SDF4"),
#'                                 tissueSiteDetail="Artery - Coronary",
#'                                 pathologyNotesCategories=TRUE)
#'
#' # Download with versioned and unversioned gencode Id.
#' expProfiles <- xQTLdownload_exp(c("ENSG00000141510.16","ENSG00000008130.15","ENSG00000078808"),
#'                                tissueSiteDetail="Artery - Coronary")
#' }
xQTLdownload_exp <- function(genes="", geneType="auto", tissueSiteDetail="Liver", toSummarizedExperiment=FALSE, recordPerChunk=80, pathologyNotesCategories=FALSE  ){
  gencodeId <- chromosome <- cutF <- genesUpper <- geneSymbol <- entrezGeneId <- tss <- description <- NULL
  .<-NULL
  cutNum <- recordPerChunk
  datasetId <- "gtex_v8"
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
    # gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV26)
    gencodeVersion <- "v26"
    genomeBuild <- "GRCh38/hg38"
  }else if(datasetId == "gtex_v7"){
    # data(tissueSiteDetailGTExv7)
    # data(gencodeGeneInfoV19)
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
    # gencodeGeneInfo <- data.table::copy(gencodeGeneInfoV19)
    gencodeVersion <- "v19"
    genomeBuild <- "GRCh37/hg19"
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be null!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }
  # convert tissueSiteDetail to tissueSiteDetailId:
  tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId


  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(genes, function(g){ stringr::str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else if( length(genes)==1 ){
      if( genes %in% gencodeGenetype$V26 | genes %in% gencodeGenetype$V19 ){
        geneType <- "geneCategory"
      }else{
        geneType <- "geneSymbol"
      }
    }else{
      geneType <- "geneSymbol"
    }
  }


  ############ convert genes. parameter check is unnecessary for this, because xQTLquery_gene check it internally.
  message("== Check the gene name :")
  geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType)
  # Only keep genes with non-na gencode ID
  geneInfo <- geneInfo[!is.na(gencodeId)]
  # for gene of the same name but with different gencodeID, like SHOX, ENSG00000185960.13 the name is in X, ENSG00000185960.13_PAR_Y is the name in Y. retain the gencode In x.
  geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & stringr::str_detect(chromosome, "Y"))]
  # for gene of the same name but with different gencodeID, like BTBD8 with entrezGeneId is NA, remove.
  geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & is.na(entrezGeneId))]
  # retain source is Source:NCBI
  # geneInfo[(genes %in% geneInfo[duplicated(geneSymbol), ]$genes)]

  # duplicates, only warning first duplicate:
  # test: genes = c("LYNX1", "TP53")
  # dupGeneSymbol <- unique(geneInfo[,.(geneSymbol, gencodeId)][,.(.SD[duplicated(geneSymbol)])]$geneSymbol[1])
  # if(!is.na(dupGeneSymbol)){
  #   message(" ")
  #   message("  == Gene [",dupGeneSymbol,"] has multiple gencode IDs: [",paste0(geneInfo[geneSymbol == dupGeneSymbol]$gencodeId, collapse = ", "),"]")
  #   stop("== Please remove duplicated genes, or take the unique gencode IDs as the input.")
  # }
  dupGencodeId <- unique(geneInfo[,.(geneSymbol, gencodeId)][,.(.SD[duplicated(gencodeId)])]$gencodeId[1])
  if(!is.na(dupGencodeId)){
    stop("duplicated gencodeId detected: ", paste0(dupGencodeId, collapse = "; "))
  }

  #
  if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
    stop("gene information is null.")
  }
  message("== Done.")
  # geneInfo <-  xQTLquery_gene(c("tp53","naDK","SDF4"), "geneSymbol", "v26", "GRCh38/hg38")

  ############ get sample info:
  message("== Get the samples' detail:")
  sampleInfo <- xQTLquery_sampleByTissue(tissueSiteDetail=tissueSiteDetail, dataType="RNASEQ", pathologyNotesCategories=pathologyNotesCategories )
  message("== Done.")
  if( !exists("sampleInfo") ||is.null(sampleInfo) ){
    stop("Failed to fetch sample information.")
  }

  ############### 分批下载：
  message("== Get the gene expression:")
  geneInfo <- cbind(geneInfo, data.table::data.table(cutF = cut(1:nrow(geneInfo),breaks=seq(0,nrow(geneInfo)+cutNum,cutNum) ) ))
  data.table::setDT(geneInfo)
  geneInfo$genesUpper <- toupper(geneInfo$genes)
  genesURL <- geneInfo[,.(genesURL=paste0(gencodeId[!is.na(gencodeId)],collapse = ",")),by=c("cutF")]
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
  for(i in 1:nrow(genesURL)){
    # construct url:
    url1 <- paste0("https://gtexportal.org/api/v2/expression/geneExpression?",
                   "datasetId=", datasetId,"&",
                   "gencodeId=", stringr::str_replace_all(genesURL[i,]$genesURL, ",","&gencodeId=" ), "&",
                   "tissueSiteDetailId=", tissueSiteDetailId
    )
    url1 <- utils::URLencode(url1)
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")

    tmp <- data.table::as.data.table(url1GetText2Json$data)
    if(nrow(tmp)==0){
      message("No expression profiles were found in ",tissueSiteDetail, " of thess ", length(genes), " genes!")
      message("== Done.")
      return(data.table::data.table())
    }
    tmp$genesUpper <- toupper(unlist(tmp[,geneType,with=FALSE]))
    # keep raw genes:
    tmp <- merge(geneInfo[cutF==genesURL[i,]$cutF,.(genes, gencodeId, geneSymbol, datasetId)], tmp[,-c("geneSymbol", "datasetId")], by="gencodeId",all.x=TRUE, sort = FALSE)
    outInfo <- rbind(outInfo, tmp)
    message("Downloaded part ", i, "/",nrow(genesURL),"; ", nrow(outInfo), " records.")
    rm(url1, url1GetText2Json, tmp)
  }
  # expression, fill with NA if NULL:
  tmpExp <- outInfo$data
  # tmpExp <- lapply(tmpExp, function(x){ if(is.null(x)){return( data.table::data.table( matrix(rep(NA,nrow(sampleInfo)),nrow = 1) ) )}else{return( data.table::data.table( matrix(x, nrow=1) ) )} })
  # expDT <- data.table::rbindlist(tmpExp)
  tmpExp <- lapply(tmpExp, function(x){  if(is.null(x)){ rep(NA, nrow(sampleInfo) )}else{ x } })
  expDT <- data.table::as.data.table(t(data.table::as.data.table(tmpExp)))
  data.table::setDF(expDT)
  colnames(expDT) <- sampleInfo$sampleId
  rownames(expDT) <- outInfo$genes
  rm(tmpExp)
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
}



#' @title Download summary statistics data for eQTLs with a specified gene, variant, tissue or study
#' @description
#'  source of all eQTL associations is EBI eQTL category.
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueLabel (character) all supported tissues can be listed using "ebi_study_tissues"
#' @param recordPerChunk (integer) number of records fetched per request (default: 1000).
#' @param study (character) name of studies can be listed using "ebi_study_tissues". If the study is null, use all studies (Default).
#' @param withB37VariantId a logical value indicating whether to return the genome location(GTEx v7) of variants. Default: FALSE.
#' @param data_source "eQTL_catalogue"(default) or "liLab"
#' @param API_version "v1"(default) or "v2"(not working)
#' @import data.table
#' @import stringr
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' # Download all eQTL associations of MLH1-rs13315355 pair in all tissues from all studies:
#' eqtlAsso <- xQTLdownload_eqtlAllAsso(gene="MLH1", variantName = "rs13315355")
#'
#' # Download eQTL associations of gene ATP11B in CD4+ T cell from all supported studies:
#' geneAsso <- xQTLdownload_eqtlAllAsso(gene="MMP7",tissueLabel = "CD4+ T cell")
#'
#' # Download eQTL associations of gene ATP11B in Muscle - Skeletal from GTEx_V8:
#' geneAsso <- xQTLdownload_eqtlAllAsso("ATP11B", tissueLabel="Muscle - Skeletal")
#'
#' # Download eQTL associations of SNP rs11568818 in Muscle - Skeletal from GTEx_V8:
#' varAsso <- xQTLdownload_eqtlAllAsso(variantName="chr11_102530930_T_C_b38",
#'                                     tissueLabel="Muscle - Skeletal")
#'
#' # Download all eQTL associations of SNP rs11568818 in all tissues from all supported studies.
#' varAsso <- xQTLdownload_eqtlAllAsso(variantName="rs11568818")
#' }
xQTLdownload_eqtlAllAsso <- function(gene="", geneType="auto", variantName="", variantType="auto", tissueLabel="", study="", recordPerChunk=1000, withB37VariantId=FALSE, data_source="eQTL_catalogue", API_version="v1" ){
  . <- geneInfoV19 <- pos <- ref<- alt<- tissue<-NULL
  chrom <- tissue_label <- study_accession <- variantId <- variant <- gencodeId <- genes<- entrezGeneId <- chromosome<- geneSymbol<- b37VariantId <- snpId <- NULL
  tissueSiteDetail <- tissueSiteDetailId <- study_label <- sample_group<- NULL
  # gene="CYP2W1"
  # geneType="geneSymbol"
  # tissueSiteDetail="Lung"

  if(data_source=="liLab"){
    tissueDT <- tissueSiteDetailGTExv8[tissueSiteDetail== tissueLabel | tissueSiteDetailId==tissueLabel][1,]
    if(nrow(tissueDT)==0){ stop("tissue not found...")}

    # check geneType
    if( !(geneType %in% c("auto","geneSymbol", "gencodeId")) ){
      stop("Parameter \"geneType\" should be choosen from \"auto\", \"geneSymbol\", and \"gencodeId\".")
    }
    if( length(gene)==1 && gene!=""){
      # Automatically determine the type of variable:
      if(geneType=="auto"){
        if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
          geneType <- "gencodeId"
        }else{
          geneType <- "geneSymbol"
        }
      }
    }
    # merge with genes:
    message("==> ",geneType)
    if(geneType!="gencodeId"){
      message("==> Querying genes...")
      geneDT <- xQTLquery_gene(gene)
      geneDT <- geneDT[,c("gencodeId")]
    }else if(geneType=="gencodeId" & all(unlist(lapply(gene, function(g){ str_detect(g, stringr::fixed(".")) }))) ){
      geneDT <- data.table(gencodeId= gene)
    }else{
      message("==> require gencodeId...")
      geneDT <- xQTLquery_gene(gene)
      geneDT <- geneDT[,c("gencodeId")]
    }
    if(nrow(geneDT)==0){ message(0);return(NULL)}
    geneDT$url1 <- paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/eQTL_b38_rsid_gene/", tissueDT$tissueSiteDetailId, "/",geneDT$gencodeId)
    # message("==> download url: ",geneDT$url1)
    geneDT$tmpFilePath <- paste(unlist(lapply(1:nrow(geneDT), function(x){tempfile(pattern = "eqtl_")})),"eQTL_",1:nrow(geneDT),".txt", sep="")

    # download sQTL by clu:
    message("== Start downloading eQTLs...")
    for(i in 1:nrow(geneDT)){
      cat("== Downloading eQTL of ", geneDT[i,]$gencodeId,  "...")
      df <- try(suppressWarnings(utils::download.file(url = geneDT[i,]$url1,
                                                      destfile=geneDT[i,]$tmpFilePath,
                                                      quiet = TRUE )), silent=TRUE)
      if(!file.exists(geneDT[i,]$tmpFilePath)){ cat("    > Failed...") }
      cat("   > Success!")
      message("\n")
    }

    message("== combine results...")
    eQTL_summary <- data.table()
    for(i in 1:nrow(geneDT)){
      if(file.exists(geneDT[i,]$tmpFilePath)){
        eQTL_summary_i <- fread(geneDT[i,]$tmpFilePath, header=TRUE)
        if(nrow(eQTL_summary_i)>0){
          # cat(i,"-")
          eQTL_summary_i$gencodeId <- geneDT[i,]$gencodeId
          eQTL_summary <- rbind(eQTL_summary, eQTL_summary_i)
          # file.remove(geneDT[i,]$tmpFilePath)
        }
      }
    }
    if(nrow(eQTL_summary)>0){
      return(eQTL_summary)
    }else{
      return(NULL)
    }
  }

  if(data_source=="eQTL_catalogue" & API_version=="v2"){
    # https://www.ebi.ac.uk/eqtl/api/v2/datasets/QTD000118/associations?size=10&start=1&variant=chr1_109274570_A_G&gene_id=ENSG00000134243
    ebi_DT <- copy(ebi_datasets)

    # check study:
    if(study==""){
      stop("study must be specified, please select form `xQTLbiolinks::ebi_datasets`")
    }else if(length(study) ==1 && study!=""){
      if(toupper(study) %in% toupper(unique(ebi_DT$study_label)) ){
        study <- unique(ebi_DT$study_label)[ toupper(unique(ebi_DT$study_label)) == toupper(study) ]
        message("== Study [", study, "] detected...")
      }else{
        message("ID\tstudy\ttissueLabel")
        a <- unique(ebi_DT[,.(study_label , sample_group)])
        for(i in 1:nrow(a)){ message(i,"\t", paste(a[i ,.(study_label , sample_group)], collapse = " \t ")) }
        stop("== Study [",study,"] can not be correctly matched, please choose from above list: ")
      }
    }

    # check tissue:
    if(tissueLabel==""){
      stop("tissueLabel must be specified, please select form `xQTLbiolinks::ebi_datasets`")
    }else if( length(tissueLabel)==1 && tissueLabel!="" ){
      if( toupper(tissueLabel) %in% toupper(unique(ebi_DT$tissue_label)) ){
        message("== Tissue label [", tissueLabel, "] detected...")
        tissueLabel <- unique(ebi_DT$tissue_label)[ toupper(unique(ebi_DT$tissue_label)) == toupper(tissueLabel) ]
      }else{
        message("ID\tstudy\ttissueLabel")
        for(i in 1:nrow(ebi_DT)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
        stop("== tissueLabel [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
      }
    }
  }

  if(data_source=="eQTL_catalogue" & API_version=="v1"){
    ebi_ST <-copy(ebi_study_tissues)

    if( length(study) ==1 && study!="" ){
      if(toupper(study) %in% toupper(unique(ebi_ST$study_accession))){
        study <- unique(ebi_ST$study_accession)[ toupper(unique(ebi_ST$study_accession)) == toupper(study) ]
        message("== Study [", study, "] detected...")
      }else{
        message("ID\tstudy\ttissueLabel")
        for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
        stop("== Study [",study,"] can not be correctly matched, please choose from above list: ")
      }
    }

    # check tissue:
    if( length(tissueLabel)==1 && tissueLabel!="" ){
      if( toupper(tissueLabel) %in% toupper(unique(ebi_ST$tissue_label)) ){
        message("== Tissue label [", tissueLabel, "] detected...")
        tissueLabel <- unique(ebi_ST$tissue_label)[ toupper(unique(ebi_ST$tissue_label)) == toupper(tissueLabel) ]
      }else{
        message("ID\tstudy\ttissueLabel")
        for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
        stop("== tissueLabel [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
      }
    }

    # check study-tissue:
    if( length(study) ==1 && length(tissueLabel)==1 && study!="" && tissueLabel!=""){
      if(nrow( ebi_ST[study_accession == study & tissue_label==tissueLabel])==1){
        message("== Study [",study,"] -- Tissue label [",tissueLabel,"] correctly mapped..")
      }else{
        message("ID\tstudy\ttissueLabel")
        for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
        stop("== Study [",study,"] -- Tissue label [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
      }
    }



    # check geneType
    if( !(geneType %in% c("auto","geneSymbol", "gencodeId")) ){
      stop("Parameter \"geneType\" should be choosen from \"auto\", \"geneSymbol\", and \"gencodeId\".")
    }
    if( length(gene)==1 && gene!=""){
      # Automatically determine the type of variable:
      if(geneType=="auto"){
        if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
          geneType <- "gencodeId"
        }else{
          geneType <- "geneSymbol"
        }
      }
    }

    # check variantType:
    if( !(variantType %in% c("auto","snpId", "variantId")) ){
      stop("Parameter \"geneType\" should be choosen from \"auto\", \"snpId\", and \"variantId\".")
    }
    if(length(variantName)==1 && variantName!=""){
      # auto pick variantType
      if(variantType=="auto"){
        if(stringr::str_detect(variantName, stringr::regex("^rs"))){
          variantType <- "snpId"
        }else if(stringr::str_count(variantName,"_")>=3){
          variantType <- "variantId"
        }else{
          stop("Note: \"variantName\" only support dbSNP id that start with \"rs\", like: rs12596338, or variant ID like: \"chr16_57156226_C_T_b38\", \"16_57190138_C_T_b37\" ")
        }
      }
    }
    # modify variant name of variantID:
    if(length(variantName)==1 && variantName!="" && variantType == "variantId"){
      variantName <- stringr::str_remove(variantName, "_b38")
    }

    # Excessive empty variables:
    if(gene=="" && variantName=="" && tissueLabel==""){
      warning("== All associations from Study [", study,"] will be returned, this will take several days...")
      warning("== Please make sure you have enough memory...")
      message("== Specified gene, variant or tissueLabel is recommended.")
    }
    if(gene=="" && variantName=="" && study==""){
      warning("== All associations in Tissue [", tissueLabel,"] will be returned, this will take several days...")
      warning("== Please make sure you have enough memory...")
      message("== Specified gene, variant or study is recommended.")
    }
    if(gene=="" && variantName=="" && tissueLabel=="" && study==""){
      stop("== All required fileds are empty, please specify gene, variant study or tissueLabel")
    }

    # if gene is not null, check geneName and add gene unversioned ensemble name.
    if(gene !=""){
      message("== Check the gene name:")
      geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType)
      # 对于只需要一个基因的：
      geneInfo <- geneInfo[!is.na(gencodeId)]
      geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & stringr::str_detect(chromosome, "Y"))]
      geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & is.na(entrezGeneId))]
      # geneInfoV19 <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v19")
      if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
        stop("Invalid gene name or type, please correct your input, or set gene with gencodeId.")
      }else{
        geneInfo$gencodeIdUnv <-stringr::str_split(geneInfo$gencodeId, stringr::fixed("."))[[1]][1]
        message("== Done.")
      }
    }

    ##################### fetch geneInfo:
    if( variantName!="" && variantType!= "auto"){
      # construct url:
      # aa<-"https://www.ebi.ac.uk/eqtl/api/associations/rs2302765?study=GTEx_V8&gene_id=ENSG00000238917&tissue=UBER_0001211"
      url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/associations/",variantName,"?links=False",
                     ifelse(gene !="", paste0("&gene_id=",geneInfo$gencodeIdUnv),""),
                     ifelse(tissueLabel!="", paste0("&tissue=",ebi_ST[tissue_label==tissueLabel]$tissue[1]),""),
                     ifelse(study!="", paste0("&study=",study),"")
      )
    }else if(gene!="" && geneType!="auto"){
      # aa <-“https://www.ebi.ac.uk/eqtl/api/genes/ENSG00000282458/associations?variant_id=chr19_80901_G_T&tissue=UBER_0001954&study=GTEx_V8”
      url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/genes/",geneInfo$gencodeIdUnv,"/associations?links=False",
                     ifelse(tissueLabel!="", paste0("&tissue=",ebi_ST[tissue_label==tissueLabel]$tissue[1]),""),
                     ifelse(study!="", paste0("&study=",study),"")
                     # ifelse(variantName!="",paste0("&variant_id=",variantName),"")
      )
    }else if( tissueLabel!=""  ){
      # aa <- "https://www.ebi.ac.uk/eqtl/api/tissues/CL_0000235/studies/Alasoo_2018/associations"
      url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/tissues/",ebi_ST[tissue_label==tissueLabel]$tissue[1],
                     ifelse(study!="", paste0("/studies/",study),""),
                     "/associations?links=False"
      )
    }else if( study!="" ){
      # aa <- "https://www.ebi.ac.uk/eqtl/api/studies/GTEx_V8/associations"
      url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/studies/", study, "/associations?links=False")
    }

    # for brain tissues with duplicated tissue id, "Brain - Cerebellar Hemisphere" and "Brain - Cerebellum"
    if( tissueLabel == "Brain - Cerebellar Hemisphere" ){
      qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellar_Hemisphere")
    }else if(tissueLabel == "Brain - Cerebellum"){
      qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellum")
    }else{
      qtl_groupStr <- ""
    }
    url1 <- paste0(url1, qtl_groupStr)
    # message(url1)

    message("== Start fetching associations...")
    # construct url:

    # # check network:
    # bestFetchMethod <- apiEbi_ping()
    # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    #   message("Note: EBI API server is busy or your network has latency, please try again later.")
    #   return(NULL)
    # }
    # gtexAsoo <- fetchContentEbi(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2],  termSize=1000)
    gtexAsoo <- fetchContentEbi(url1, method = "fromJSON", downloadMethod = "auto",  termSize=recordPerChunk)
    if(is.null(gtexAsoo)){
      return(NULL)
    }
    gtexAsooList <- do.call(c, gtexAsoo)

    if(length(gtexAsooList)==0){
      message("No association found!")
      return(NULL)
    }
    gtexAsooList <- lapply(gtexAsooList, function(x){ data.table(variantId= x$variant, snpId=x$rsid,type=x$type,maf=x$maf,beta=x$beta,
                                                                 chrom=x$chromosome, pos=x$position,ref=x$ref, alt=x$alt,study_id=x$study_id,
                                                                 se=x$se, median_tpm=x$median_tpm, pValue=x$pvalue, totalAlleles=x$an, allelCounts=x$ac, imputationR2=x$r2,
                                                                 tissue_label=x$tissue_label,tissue=x$tissue,qtl_group=x$qtl_group, condition=x$condition,
                                                                 molecular_trait_id = x$molecular_trait_id, gene_id= x$gene_id
    ) })
    gtexAsooDT <- data.table::rbindlist(gtexAsooList, fill=TRUE)
    gtexAsooDT[,variantId:=.(paste0( "chr",stringr::str_remove_all(chrom, "chr"),"_",pos,"_",ref,"_",alt ))]
    if( !("tissue_label" %in% names(gtexAsooDT)) ){
      if(tissueLabel==""){
        ebi_ST_tmp <- ebi_ST[tissue_label!="Brain - Cerebellum", .(tissue, tissue_label,study_id=study_accession)]
      }else{ ebi_ST_tmp <-ebi_ST[tissue_label==tissueLabel,.(tissue, tissue_label,study_id=study_accession)]  }
      gtexAsooDT <- merge(gtexAsooDT, ebi_ST_tmp, by = c("tissue","study_id"))
    }
    gtexAsooDT <- cbind(gtexAsooDT[,-c("study_id", "tissue", "tissue_label", "qtl_group")], gtexAsooDT[,c("study_id", "tissue", "tissue_label", "qtl_group")])
    if(gene!="" && geneType!="auto"){
      gtexAsooDT$geneSymbol <- geneInfo$geneSymbol
      gtexAsooDT$gencodeId_GTEX_v8 <- geneInfo$gencodeId
      # gtexAsooDT$gencodeId_GTEX_v7 <- ifelse( exists("geneInfoV19") && nrow(geneInfoV19)>0, geneInfoV19$gencodeId, "")
    }else{
      geneInfo <- xQTLquery_gene(unique(gtexAsooDT$gene_id))
      gtexAsooDT <- merge(gtexAsooDT[,], geneInfo[,.(gene_id=genes, geneSymbol, gencodeId_GTEX_v8=gencodeId)], by = "gene_id", all.x = TRUE)[,-c("gene_id")]
    }
    if(withB37VariantId){
      # add dbSNP id and  hg19 cordinate:
      gtexAsooDTb37 <- xQTLquery_varPos(chrom = paste0("chr",unique(gtexAsooDT$chrom)), pos = unique(gtexAsooDT$pos), recordPerChunk = 250)
      gtexAsooDTb37$variantId <- unlist(lapply(gtexAsooDTb37$variantId, function(x){ splitInfo=stringr::str_split(x, stringr::fixed("_"))[[1]]; paste0(splitInfo[-5], collapse="_") }))
      gtexAsooDT <- merge(gtexAsooDT, gtexAsooDTb37[,.(variantId, b37VariantId)], by=c("variantId"), all.x=TRUE )
      # gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
      # gtexAsooDT <- cbind(gtexAsooDT[,.(snpId, variantId, b37VariantId)], gtexAsooDT[,-c("snpId", "variantId", "b37VariantId", "chrom", "pos")])
    }
    gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
    return(gtexAsooDT)
  }
}



#' @title Download summary statistics data for eQTLs with genome positions.
#'
#' @param chrom (character) name of chromesome, including chr1-chr22, chrX.
#' @param pos_lower (integer) lower base pair location threshold, expressed as an integer
#' @param pos_upper (integer) upper base pair location threshold, expressed as an integer
#' @param p_lower (numeric) lower p-value threshold, can be expressed as a float or using mantissa and exponent annotation (0.001 or 1e-3 or 1E-3)
#' @param p_upper (numeric) upper p-value threshold, can be expressed as a float or using mantissa and exponent annotation (0.001 or 1e-3 or 1E-3)
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueLabel (character) all supported tissues can be listed using "ebi_study_tissues".
#' @param study (character) name of studies can be listed using "ebi_study_tissues".
#' @param recordPerChunk (integer) number of records fetched per request (default: 1000).
#' @param withB37VariantId a logical value indicating whether to return the genome location(GTEx v7) of variants. Default: FALSE.
#'
#' @import data.table
#' @import stringr
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#' eqtlAssos <- xQTLdownload_eqtlAllAssoPos(chrom = "chr11",
#'                                         pos_lower=101400000, pos_upper = 101400013,
#'                                         tissueLabel="Brain - Cerebellar Hemisphere",
#'                                         )
#' }
xQTLdownload_eqtlAllAssoPos <- function(chrom="", pos_lower=numeric(0), pos_upper=numeric(0), p_lower=0, p_upper=1.1,  gene="", geneType="auto", tissueLabel="", study="gtex_v8", recordPerChunk=1000, withB37VariantId=FALSE){
  .<-NULL
  study_accession <- tissue_label <- gencodeId <- genes <- geneSymbol <- chromosome <- entrezGeneId <- variantId <- pos <- ref <- alt <- tissue <- b37VariantId <- NULL
  ebi_ST <- data.table::copy(ebi_study_tissues)

  # check chrome:
  if(length(chrom)==1 && chrom!="" && !is.na(chrom) ){
    chrom <- toupper(stringr::str_remove(chrom, "chr"))
  }else if( !(chrom %in% c(1:22, "X")) ){
    stop("Chromosome must be choosen from one of: ", paste0(paste0("chr",c(1:22, "X")), collapse = ", "))
  }else{
    message("== chromosome: ", chrom)
  }

  # check position and pvalue:
  if( length(pos_lower)==1 && length(pos_upper) ==1 && length(p_upper)==1 && length(p_lower)==1 && is.wholenumber(pos_lower) && is.wholenumber(pos_upper) && is.numeric(p_lower) && is.numeric(p_upper)){
    message("== Genome range: ",chrom,":", as.character(as.integer(pos_lower)),"-",as.character(as.integer(pos_upper)))
  }else{
    stop("pos_lower, pos_upper must be the whole numbers, and p_lower, p_upper must be the numeric value.")
  }
  if(pos_upper < pos_lower){stop("pos_upper must greater than pos_lower")}
  if( p_upper < p_lower ){ stop("p_upper must greater than p_lower")}


  # check study:
  if( length(study) ==1 && study!="" ){
    if(toupper(study) %in% toupper(unique(ebi_ST$study_accession))){
      study <- unique(ebi_ST$study_accession)[ toupper(unique(ebi_ST$study_accession)) == toupper(study) ]
      message("== Study [", study, "] detected...")
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== Study [",study,"] can not be correctly matched, please choose from above list: ")
    }
  }

  # check tissue:
  if( length(tissueLabel)==1 && tissueLabel!="" ){
    if( toupper(tissueLabel) %in% toupper(unique(ebi_ST$tissue_label)) ){
      message("== Tissue label [", tissueLabel, "] detected...")
      tissueLabel <- unique(ebi_ST$tissue_label)[ toupper(unique(ebi_ST$tissue_label)) == toupper(tissueLabel) ]
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== tissueLabel [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
    }
  }

  # check study-tissue:
  if( length(study) ==1 && length(tissueLabel)==1 && study!="" && tissueLabel!=""){
    if(nrow( ebi_ST[study_accession == study & tissue_label==tissueLabel])==1){
      message("== Study [",study,"] -- Tissue label [",tissueLabel,"] correctly mapped..")
    }else{
      message("ID\tstudy\ttissueLabel")
      for(i in 1:nrow(ebi_ST)){ message(i,"\t", paste(ebi_study_tissues[i ,.(study_accession, tissue_label)], collapse = " \t ")) }
      stop("== Study [",study,"] -- Tissue label [",tissueLabel,"] can not be correctly matched, please choose from above list: ")
    }
  }

  # check geneType
  if( !(geneType %in% c("auto","geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"auto\", \"geneSymbol\", and \"gencodeId\".")
  }
  if( length(gene)==1 && gene!=""){
    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }
  }

  # if gene is not null, check geneName and add gene unversioned ensemble name.
  if(gene !=""){
    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType)
    # 对于只需要一个基因的：
    geneInfo <- geneInfo[!is.na(gencodeId)]
    geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & stringr::str_detect(chromosome, "Y"))]
    geneInfo <- geneInfo[!(genes %in% geneInfo[duplicated(geneSymbol), ]$genes & is.na(entrezGeneId))]
    # geneInfoV19 <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = "v19")
    if(nrow(geneInfo)==0 || is.null(geneInfo)|| !exists("geneInfo") ){
      stop("Invalid gene name or type, please correct your input, or set gene with gencodeId.")
    }else{
      geneInfo$gencodeIdUnv <-stringr::str_split(geneInfo$gencodeId, stringr::fixed("."))[[1]][1]
      message("== Done.")
    }
  }

  message("== Start fetching associations...", format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))
  # url1 <- "https://www.ebi.ac.uk/eqtl/api/chromosomes/11/associations?study=GTEx_V8&gene_id=ENSG00000137673&bp_lower=101798614&bp_upper=103462313&p_upper=1e-1"
  # url1 <- "https://www.ebi.ac.uk/eqtl/api/v2/datasets/QTD000266/associations?pos=1:123456-124456"

  url1 <- paste0("https://www.ebi.ac.uk/eqtl/api/chromosomes/", chrom, "/associations?links=False",
                 "&bp_lower=",as.character(as.integer(pos_lower)), "&bp_upper=",as.character(as.integer(pos_upper)), "&p_upper=",p_upper, "&p_lower=",p_lower,
                 ifelse(gene !="", paste0("&gene_id=",geneInfo$gencodeIdUnv),""),
                 ifelse(tissueLabel!="", paste0("&tissue=",ebi_ST[tissue_label==tissueLabel]$tissue[1]),""),
                 ifelse(study!="", paste0("&study=",study),"")
                 )
  message(url1)
  # for brain tissues with duplicated tissue id, "Brain - Cerebellar Hemisphere" and "Brain - Cerebellum"
  if( tissueLabel == "Brain - Cerebellar Hemisphere" ){
    qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellar_Hemisphere")
  }else if(tissueLabel == "Brain - Cerebellum"){
    qtl_groupStr <- paste0("&qtl_group=Brain_Cerebellum")
  }else{
    qtl_groupStr <- ""
  }
  url1 <- paste0(url1, qtl_groupStr)

  # message(url1)
  gtexAsoo <- fetchContentEbi(url1, method = "fromJSON", downloadMethod = "auto",  termSize=recordPerChunk)
  gtexAsooList <- do.call(c, gtexAsoo)
  if(length(gtexAsooList)==0){
    message("No association found!")
    return(NULL)
  }
  gtexAsooList <- lapply(gtexAsooList, function(x){ data.table(variantId= x$variant, snpId=x$rsid,type=x$type,maf=x$maf,beta=x$beta,
                                                               chrom=x$chromosome, pos=x$position,ref=x$ref, alt=x$alt,study_id=x$study_id,
                                                               se=x$se, median_tpm=x$median_tpm, pValue=x$pvalue, totalAlleles=x$an, allelCounts=x$ac, imputationR2=x$r2,
                                                               tissue_label=x$tissue_label,tissue=x$tissue,qtl_group=x$qtl_group, condition=x$condition,
                                                               molecular_trait_id = x$molecular_trait_id, gene_id= x$gene_id
  ) })
  gtexAsooDT <- data.table::rbindlist(gtexAsooList, fill=TRUE)
  gtexAsooDT[,variantId:=.(paste0( "chr",stringr::str_remove_all(chrom, "chr"),"_",pos,"_",ref,"_",alt ))]
  if( !("tissue_label" %in% names(gtexAsooDT)) ){
    if(tissueLabel==""){
      ebi_ST_tmp <- ebi_ST[tissue_label!="Brain - Cerebellum", .(tissue, tissue_label,study_id=study_accession)]
    }else{ ebi_ST_tmp <-ebi_ST[tissue_label==tissueLabel,.(tissue, tissue_label,study_id=study_accession)]  }
    gtexAsooDT <- merge(gtexAsooDT, ebi_ST_tmp, by = c("tissue","study_id"))
  }
  gtexAsooDT <- cbind(gtexAsooDT[,-c("study_id", "tissue", "tissue_label", "qtl_group")], gtexAsooDT[,c("study_id", "tissue", "tissue_label", "qtl_group")])

  # add gene info.
  if(gene!="" && geneType!="auto"){
    gtexAsooDT$geneSymbol <- geneInfo$geneSymbol
    gtexAsooDT$gencodeId_GTEX_v8 <- geneInfo$gencodeId
    # gtexAsooDT$gencodeId_GTEX_v7 <- ifelse( exists("geneInfoV19") && nrow(geneInfoV19)>0, geneInfoV19$gencodeId, "")
  }else{
    geneInfo <- xQTLquery_gene(unique(gtexAsooDT$gene_id))
    gtexAsooDT <- merge(gtexAsooDT[,], geneInfo[,.(gene_id=genes, geneSymbol, gencodeId_GTEX_v8=gencodeId)], by = "gene_id", all.x = TRUE)[,-c("gene_id")]
  }
  if(withB37VariantId){
    # add dbSNP id and  hg19 cordinate:
    gtexAsooDTb37 <- xQTLquery_varPos(chrom = paste0("chr",unique(gtexAsooDT$chrom)), pos = unique(gtexAsooDT$pos), recordPerChunk = 250)
    gtexAsooDTb37$variantId <- unlist(lapply(gtexAsooDTb37$variantId, function(x){ splitInfo=stringr::str_split(x, stringr::fixed("_"))[[1]]; paste0(splitInfo[-5], collapse="_") }))
    gtexAsooDT <- merge(gtexAsooDT, gtexAsooDTb37[,.(variantId, b37VariantId)], by=c("variantId"), all.x=TRUE )
    # gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
    # gtexAsooDT <- cbind(gtexAsooDT[,.(snpId, variantId, b37VariantId)], gtexAsooDT[,-c("snpId", "variantId", "b37VariantId", "chrom", "pos")])
  }
  gtexAsooDT$variantId <- paste0(gtexAsooDT$variantId,"_b38")
  return(gtexAsooDT)
}


#' @title Download normalized expression for gene with a variant-gene pair
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @import data.table
#' @import stringr
#' @import jsonlite
#' @import stats
#' @return A data.table object.
#' @export
#' @examples
#' \donttest{
#' # Download exp with variant-gene pair in different tissues:
#' xQTLdownload_eqtlExp(variantName="rs1641513",gene="TP53", tissueSiteDetail="Liver")
#'
#' # Download expression using variant ID and gencode ID.
#' xQTLdownload_eqtlExp(variantName="chr1_14677_G_A_b38",gene="ENSG00000228463.9",
#'                      tissueSiteDetail="Stomach")
#' }
xQTLdownload_eqtlExp <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail=""){
  gencodeGenetype <- chromosome <-gencodeId <-NULL
  .<-NULL
  datasetId="gtex_v8"

  # variantName="chr1_14677_G_A_b38"
  # variantType="variantId"
  # datasetId="gtex_v8"
  # tissueSiteDetail="Liver"
  # gene = "ENSG00000228463.9"
  # geneType="gencodeId"
  # check version:
  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if( datasetId == "gtex_v7" ){
    gencodeVersion <- "v19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }
  # check gene:
  if( length(gene) >1 || !is.character(gene) ){
    stop("Parameter \"gene\" must be a character string. ")
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if(  !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }else{
    message("Invalid parameter \"tissueSiteDetail\".")
    return(data.table::data.table())
  }

  ##################### fetch varInfo
  if(variantName!=""){

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

    message("== Check the variant name:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType)
    if(nrow(varInfo)==0 || is.null(varInfo) || !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"variantName\" can not be null!")
  }

  ##################### fetch geneInfo:
  if(gene !=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, recordPerChunk = 150)
    geneInfo <- geneInfo[!is.na(gencodeId)]
    if(nrow(geneInfo)==0|| is.null(geneInfo) || !exists("geneInfo")){
      stop("Invalid gene name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"gene\" can not be null!")
  }

  message("== Downloading expression...")
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # construct url
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/api/v2/association/dyneqtl?",
                 "gencodeId=",geneInfo$gencodeId,"&",
                 "variantId=",varInfo$variantId,"&",
                 "tissueSiteDetailId=",tissueSiteDetailId,"&",
                 "datasetId=",datasetId
  )
  url1 <- utils::URLencode(url1)
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  expData <- data.table::data.table(normExp=url1GetText2Json$data,genotypes=url1GetText2Json$genotypes)

  if(nrow(expData)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & gene!="","-",""),ifelse(gene=="","",paste0(" gene [", gene,"]")), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  message("=================================")
  message("== Summary: ")
  message("[ pValue ]:          ",url1GetText2Json$pValue)
  message("[ pValueThreshold ]: ",url1GetText2Json$pValueThreshold)
  message("[ nes ]:             ",url1GetText2Json$nes)
  message("[ maf ]:             ",url1GetText2Json$maf)
  message("[ error ]:           ",url1GetText2Json$error)
  message("== Normalized expression and genotypes of [", nrow(expData), "] samples were found for", paste0(" variant [", variantName,"]"), " -",paste0(" gene [", gene,"]")," pair in tissue [", tissueSiteDetail,"] in [",datasetId,"]."  )
  message("== Genotype: ",url1GetText2Json$hetCount," het; ",url1GetText2Json$homoAltCount," hom; " ,url1GetText2Json$homoRefCount," ref.")

  message("== For more normalization method details, you can visit: https://www.gtexportal.org/home/faq#normalization. ")
  return(expData)
}


#' @title Download normalized PSI value of intron for a sQTL pair
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param phenotypeId A character string. Format like: "chr1:497299:498399:clu_54863:ENSG00000239906.1"
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @import data.table
#' @import stringr
#' @import jsonlite
#' @import stats
#' @return A data.table object.
#' @export
#'
#' @examples
#' # Download sQTL expression in different tissues:
#' xQTLdownload_sqtlExp(variantName="rs1450891501",
#'                      phenotypeId="chr1:497299:498399:clu_54863:ENSG00000239906.1",
#'                      tissueSiteDetail="Lung")
#'
#' # Dowload sQTL expression using variant ID.
#' xQTLdownload_sqtlExp(variantName="chr1_14677_G_A_b38",
#'                      phenotypeId="chr1:15947:16607:clu_40980:ENSG00000227232.5",
#'                      tissueSiteDetail="Whole Blood")
xQTLdownload_sqtlExp <- function(variantName="", phenotypeId="", variantType="auto", tissueSiteDetail=""){
  # variantName="chr1_739465_TTTTG_T_b38"
  # phenotypeId="chr1:497299:498399:clu_54863:ENSG00000239906.1"
  # variantType="variantId"
  # tissueSiteDetail="Lung"
  # datasetId="gtex_v8"
  # check version:
  datasetId="gtex_v8"
  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if( datasetId == "gtex_v7" ){
    gencodeVersion <- "v19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }
  # check gene:
  if( length(phenotypeId) >1 || !is.character(phenotypeId) ){
    stop("Parameter \"phenotypeId\" must be a character string. ")
  }
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if(  !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }else{
    message("Invalid parameter \"tissueSiteDetail\".")
    return(data.table::data.table())
  }

  ##################### fetch varInfo
  if(variantName!=""){

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

    message("== Check the variant name:")
    varInfo <- xQTLquery_varId(variantName=variantName, variantType = variantType)
    if(nrow(varInfo)==0 || is.null(varInfo) || !exists("varInfo")){
      stop("Invalid variant name or type, please correct your input.")
    }else{
      message("== Done.")
    }
  }else{
    stop("Parameter \"variantName\" can not be null!")
  }

  ##################### fetch geneInfo:
  # if(gene !=""){
  #   message("== Querying gene info from API server:")
  #   geneInfo <- xQTLquery_gene(genes=gene, geneType = geneType, gencodeVersion = gencodeVersion, recordPerChunk = 150)
  #   if(nrow(geneInfo)==0|| is.null(geneInfo) || !exists("geneInfo")){
  #     stop("Invalid gene name or type, please correct your input.")
  #   }else{
  #     message("== Done.")
  #   }
  # }else{
  #   stop("Parameter \"gene\" can not be null!")
  # }

  message("== Downloading expression...")
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # construct url
  ########## construct url for sig association
  url1 <- paste0("https://gtexportal.org/api/v2/association/dynsqtl?",
                 "phenotypeId=",phenotypeId,"&",
                 "variantId=",varInfo$variantId,"&",
                 "tissueSiteDetailId=",tissueSiteDetailId,"&",
                 "datasetId=",datasetId
  )
  url1 <- utils::URLencode(url1)
  # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  expData <- data.table::data.table(normExp=url1GetText2Json$data,genotypes=url1GetText2Json$genotypes)

  if(nrow(expData)==0){
    message("No significant associations were found for", ifelse(variantName=="","",paste0(" variant [", variantName,"]")), ifelse(variantName!="" & phenotypeId!="","-",""),ifelse(phenotypeId=="","",paste0(" phenotypeId [", phenotypeId,"]")), " in ",datasetId)
    return(data.table::data.table())
  }else{
    message("== Done.")
  }
  message("=================================")
  message("== Summary: ")
  message("[ pValue ]:          ",url1GetText2Json$pValue)
  message("[ pValueThreshold ]: ",url1GetText2Json$pValueThreshold)
  message("[ nes ]:             ",url1GetText2Json$nes)
  message("[ maf ]:             ",url1GetText2Json$maf)
  message("[ error ]:           ",url1GetText2Json$error)
  message("== Normalized expression and genotypes of [", nrow(expData), "] samples were found for", paste0(" variant [", variantName,"]"), " -",paste0(" phenotypeId [", phenotypeId,"]")," pair in tissue [", tissueSiteDetail,"] in [",datasetId,"]."  )
  message("== Genotype: ",url1GetText2Json$hetCount," het; ",url1GetText2Json$homoAltCount," hom; " ,url1GetText2Json$homoRefCount," ref.")

  message("== For more normalization method details, you can visit: https://www.gtexportal.org/home/faq#normalization. ")
  return(expData)
}



#' @title Download eGenes (eQTL Genes) for a specified gene or a tissue
#' @description
#'  eGenes are genes that have at least one significant cis-eQTL acting upon them. Results can be filtered by tissue.
#' @param gene  (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param recordPerChunk (integer) number of records fetched per request (default: 2000).
#' @import data.table
#' @import stringr
#' @import utils
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' eGeneInfo <- xQTLdownload_egene(tissueSiteDetail="Prostate", gene="CICP3")
#' eGeneInfo <- xQTLdownload_egene(tissueSiteDetail="Prostate")
#' }
xQTLdownload_egene <- function(tissueSiteDetail="", gene = "", geneType="auto", recordPerChunk=2000){
  .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss <- log2AllelicFoldChange <- empiricalPValue <- pValue <- pValueThreshold <- qValue <-NULL
  # gene="DDX11"
  # geneType="geneSymbol"
  # datasetId = "gtex_v8"
  # tissueSiteDetail="Lung"
  # recordPerChunk=100
  datasetId = "gtex_v8"

  page_tmp <- 0
  pageSize_tmp <- recordPerChunk

  # parameter check: datasetId
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail))  ){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }

  # Fetch gene info:
  if(gene!=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType)
    geneInfo <- geneInfo[!is.na(gencodeId)]
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",gene,"] you entered could not be found!")
    }
    message("== Done.")
  }

  message("Downloading eGenes..")
  # url1 <- "https://gtexportal.org/api/v2/association/egene?itemsPerPage=250&searchTerm=ENSG00000013573.16&sortBy=log2AllelicFoldChange&sortDirection=asc&tissueSiteDetailId=Thyroid&datasetId=gtex_v8"
  outInfo <- data.table::data.table()
  url1 <- paste0("https://gtexportal.org/api/v2/association/egene?",
                 "page=",page_tmp,"&",
                 "itemsPerPage=",pageSize_tmp,
                 ifelse(gene=="","",paste0("&searchTerm=",geneInfo$gencodeId)),
                 "&sortBy=log2AllelicFoldChange&sortDirection=asc",
                 ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                 "&datasetId=",datasetId)
  url1 <- utils::URLencode(url1)
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")
  # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  tmp <- data.table::as.data.table(url1GetText2Json$data)
  outInfo <- rbind(outInfo, tmp)
  message("Records: ",nrow(outInfo),"/",url1GetText2Json$paging_info$totalNumberOfItems,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$paging_info$numberOfPages)
  page_tmp<-page_tmp+1
  while( page_tmp <= (url1GetText2Json$paging_info$numberOfPages-1) ){
    url1 <- paste0("https://gtexportal.org/api/v2/association/egene?",
                   "page=",page_tmp,"&",
                   "itemsPerPage=",pageSize_tmp,
                   ifelse(gene=="","&",paste0("&searchTerm=",geneInfo$gencodeId)),
                   "sortBy=log2AllelicFoldChange&sortDirection=asc",
                   ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                   "&datasetId=",datasetId)
    url1 <- utils::URLencode(url1)
    # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    tmp <- data.table::as.data.table(url1GetText2Json$data)
    outInfo <- rbind(outInfo, tmp)
    message("Records: ",nrow(outInfo),"/",url1GetText2Json$paging_info$totalNumberOfItems,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$paging_info$numberOfPages)
    page_tmp <- page_tmp+1
  }
  outInfo <- merge(outInfo, tissueSiteDetailGTEx, by="tissueSiteDetailId")

  if(gene!=""){
    outInfo <- merge( geneInfo[,.(gencodeId, geneSymbol, entrezGeneId, geneType, chromosome, start, end, tss)],
                      outInfo[,.(gencodeId, log2AllelicFoldChange, empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)],
                      by="gencodeId")
  }else{
    outInfo <- outInfo[,.(gencodeId, log2AllelicFoldChange, empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)]
  }

  return(outInfo)

}


#' @title Download details of sGenes (sQTL Genes) for a specified gene or a tissue.
#' @description
#'  sGenes are genes that have at least one significant sQTL acting upon them. Results may be filtered by tissue.
#' @param gene  (character) gene symbol or gencode id (versioned or unversioned are both supported). Can be null.
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param recordPerChunk (integer) number of records fetched per request (default: 2000).
#' @import data.table
#' @import stringr
#' @import utils
#' @return A data.table object.
#' @export
#'
#' @examples
#' \donttest{
#' sGeneInfo <- xQTLdownload_sgene(tissueSiteDetail="Liver")
#' sGeneInfo <- xQTLdownload_sgene(tissueSiteDetail="Liver", gene="DDX11" )
#' }
xQTLdownload_sgene <- function( tissueSiteDetail="", gene = "", geneType="auto", recordPerChunk=2000){
  phenotypeId <- nPhenotypes <- .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss <- log2AllelicFoldChange <- empiricalPValue <- pValue <- pValueThreshold <- qValue <-NULL
  datasetId = "gtex_v8"
  # gene="ENSG00000013573.16"
  # geneType="gencodeId"
  # datasetId = "gtex_v8"
  # tissueSiteDetail="Thyroid"
  # recordPerChunk=2000

  page_tmp <- 0
  pageSize_tmp <- recordPerChunk


  # parameter check: datasetId
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail))){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }

  # Fetch gene info:
  if(gene!=""){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Check the gene name:")
    geneInfo <- xQTLquery_gene(genes=gene, geneType=geneType)
    geneInfo <- geneInfo[!is.na(gencodeId)]
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",gene,"] you entered could not be found!")
    }
    message("== Done.")
  }

  message("Downloading sgenes...")

  # url1 <- "https://gtexportal.org/api/v2/association/egene?itemsPerPage=250&searchTerm=ENSG00000013573.16&sortBy=log2AllelicFoldChange&sortDirection=asc&tissueSiteDetailId=Thyroid&datasetId=gtex_v8"
  outInfo <- data.table::data.table()
  url1 <- paste0("https://gtexportal.org/api/v2/association/sgene?",
                 "page=",page_tmp,"&",
                 "itemsPerPage=",pageSize_tmp,
                 "&sortBy=empiricalPValue&sortDirection=asc&",
                 "datasetId=",datasetId,
                 ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                 ifelse(gene=="","",paste0("&gencodeId=",geneInfo$gencodeId))
                 )
  # message(url1)
  url1 <- utils::URLencode(url1)
  # bestFetchMethod <- apiAdmin_ping()
  # if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
  #   # message("Note: API server is busy or your network has latency, please try again later.")
  #   return(NULL)
  # }
  # message("GTEx API successfully accessed!")
  # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
  tmp <- data.table::as.data.table(url1GetText2Json$data)
  outInfo <- rbind(outInfo, tmp)
  message("Records: ",nrow(outInfo),"/",url1GetText2Json$paging_info$totalNumberOfItems,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$paging_info$numberOfPages)
  page_tmp<-page_tmp+1
  while( page_tmp <= (url1GetText2Json$paging_info$numberOfPages-1) ){
    url1 <- paste0("https://gtexportal.org/api/v2/association/sgene?",
                   "page=",page_tmp,"&",
                   "itemsPerPage=",pageSize_tmp,
                   "&sortBy=empiricalPValue&sortDirection=asc&",
                   "datasetId=",datasetId,
                   ifelse(tissueSiteDetail=="","&",paste0("&tissueSiteDetailId=",tissueSiteDetailId)),
                   ifelse(gene=="","",paste0("&gencodeId=",geneInfo$gencodeId))
    )
    url1 <- utils::URLencode(url1)
    # suppressMessages(url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2]))
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    tmp <- data.table::as.data.table(url1GetText2Json$sgene)
    outInfo <- rbind(outInfo, tmp)
    message("Records: ",nrow(outInfo),"/",url1GetText2Json$paging_info$totalNumberOfItems,"; downloaded: ", page_tmp+1, "/", url1GetText2Json$paging_info$numberOfPages)
    page_tmp <- page_tmp+1
  }
  outInfo <- merge(outInfo, tissueSiteDetailGTEx, by="tissueSiteDetailId")

  if(gene!=""){
    outInfo <- merge( geneInfo[,.(gencodeId, geneSymbol, entrezGeneId, geneType, chromosome, start, end, tss)],
                      outInfo[,.(gencodeId, phenotypeId, nPhenotypes, empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)],
                      by="gencodeId")
  }else{
    outInfo <- outInfo[,.(gencodeId, phenotypeId, nPhenotypes,empiricalPValue, pValue, pValueThreshold, qValue, tissueSiteDetail, datasetId)]
  }

  return(outInfo)

}


#' @title Download median expressions for multiple genes in a specified tissue
#' @param genes (character string or a character vector) gene symbols or gencode ids (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param recordPerChunk (integer) number of records fetched per request (default: 150).
#' @import data.table
#' @import utils
#' @import stringr
#' @return A data.table object.
#' @export
#'
#' @examples
#' geneMedExp <- xQTLdownload_geneMedExp(genes="LYNX1")
#' geneMedExp <- xQTLdownload_geneMedExp(genes=c("TP53", "IRF5"))
xQTLdownload_geneMedExp <- function(genes="", geneType="auto", tissueSiteDetail="", recordPerChunk=150 ){
  .<-NULL
  gencodeId <- geneSymbol <- entrezGeneId <- chromosome <- tss<-strand <- NULL
  datasetId="gtex_v8"
  # check genes
  if( is.null(genes) ||  any(is.na(genes)) || any(genes=="") ||length(genes)==0 ){
    stop("Parameter \"genes\" can not be NULL or NA!")
  }

  # geneType
  if( is.null(geneType) ||  any(is.na(geneType)) || any(geneType=="") || length(geneType)!=1){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(geneType %in% c("auto","geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"auto\",\"geneSymbol\", \"gencodeId\".")
  }


  # parameter check: datasetId
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail))){
    stop("Parameter \"tissueSiteDetail\" can not be NULL or NA!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be a character string!")
  }else if( tissueSiteDetail!="" && !(tissueSiteDetail %in% c( tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c("", paste0(1:nrow(tissueSiteDetailGTEx),". ",tissueSiteDetailGTEx$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }else if( tissueSiteDetail!="" ){
    # convert tissueSiteDetail to tissueSiteDetailId:
    tissueSiteDetailId <- tissueSiteDetailGTEx[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId
  }

  # Fetch gene info:
  if(all(genes!="") || length(genes)>0){

    # Automatically determine the type of variable:
    if(geneType=="auto"){
      if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
        geneType <- "gencodeId"
      }else{
        geneType <- "geneSymbol"
      }
    }

    message("== Check the gene name:")
    suppressMessages( geneInfo <- xQTLquery_gene(genes=genes, geneType=geneType, recordPerChunk = recordPerChunk))
    geneInfo <- geneInfo[!is.na(gencodeId)]
    if(nrow(geneInfo)==0 || is.null(geneInfo)||!exists("geneInfo") ){
      stop("The gene [",stringr::str_sub(paste(genes, collapse = ","),1,20),"....] you entered could not be found!")
    }
    message("== Done.")
  }

  message("Downloading gene median expression...",format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))


  #
  outInfo <- data.table::data.table()
  genesCut <- data.table::data.table(gencodeId=geneInfo$gencodeId, ID=1:nrow(geneInfo), cutF = as.character(cut(1:nrow(geneInfo),breaks=seq(0,nrow(geneInfo)+recordPerChunk,recordPerChunk) )) )
  genesURL <- genesCut[,.(genesURL=paste0(gencodeId,collapse = ",")),by=c("cutF")]
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
  for(i in 1:nrow(genesURL)){
    # construct url:
    url1 <- paste0("https://gtexportal.org/api/v2/expression/medianGeneExpression?",
                   "datasetId=",datasetId,"&",
                   "gencodeId=", str_replace_all(genesURL[i,]$genesURL,",", "&gencodeId="),
                   ifelse(tissueSiteDetail=="","&", paste0("&tissueSiteDetailId=", tissueSiteDetailId))
    )
    url1 <- utils::URLencode(url1)
    # url1GetText2Json <- fetchContent(url1, method = bestFetchMethod[1], downloadMethod = bestFetchMethod[2])
    url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto")
    url1GetText2Json2DT <- data.table::as.data.table(url1GetText2Json$data)
    if( nrow(url1GetText2Json2DT)==0 ){
      message( "0 record fatched!" )
      return(NULL)
    }
    tmp <- url1GetText2Json2DT[,.(gencodeId, geneSymbol, median, unit,tissueSiteDetailId)]

    outInfo <- rbind(outInfo, tmp)
    message("Downloaded  ", nrow(outInfo), " records.")
    # message("Downloaded  ", round(i/nrow(genesURL)*100,2),"%; totally ", length(na.omit(outInfo$gencodeId)), " records fetched!")
    # rm(url1, url1Get, url1GetText, url1GetText2Json, url1GetText2Json2DT)
  }
  outInfo <- merge(outInfo, tissueSiteDetailGTEx, by="tissueSiteDetailId")
  outInfo <- merge(outInfo[,.(gencodeId, geneSymbol, median, tissueSiteDetail)], geneInfo[,.(gencodeId, geneSymbol, entrezGeneId, geneType, chromosome, start, end, strand, tss)], by=c("gencodeId", "geneSymbol"))
  message("Unit of expression: ",unique(tmp$unit))
  return(outInfo)
}

#' @title Retrieve SNP pairwise LD from locuscompare database
#' @description
#'  SNP pairwise lD are calculated based on 1000 Genomes Project Phase 3 version 5.
#'  For storage-efficiency, the output will only include SNPs with r2 > 0.2 with the input SNP.
#' @param chr (string) Chromosome name. e.g. '22'. Notice that the name should not contain 'chr'.
#' @param snp (string) SNP rsID.
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @import RMySQL
#' @import DBI
#' @return A data.frame object.
#' @export
#' @examples
#' ld <- retrieveLD('6', 'rs9349379', 'AFR')
retrieveLD = function(chr, snp, population){
  # conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
  conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare", "locuscomparer" ,"12345678","locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com")
  on.exit(RMySQL::dbDisconnect(conn))

  chr <- str_remove_all(chr,"chr")
  res1 = DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_A, SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_A = '%s';",
      chr,
      population,
      snp
    )
  )

  res2 = DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_B as SNP_A, SNP_A as SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_B = '%s';",
      chr,
      population,
      snp
    )
  )

  res = rbind(res1,res2)
  return(res)
}



#' @title Download summary statistics data for sQTLs with a specified gene or a tissue
#'
#' @param genes (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "gencodeId".
#' @param tissue (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8`
#' @param clu_names (character) If provided, only the sQTL of clu_names will be downloaded
#' @param clu_geneid_DF (data.frame) If provided, clu-gencode mapping relationship will be loaded from this data.frame.
#' @import data.table
#' @import stringr
#' @importFrom utils download.file
#' @return A data.table object of sQTL dataset.
#' @export
#'
#' @examples
#' \donttest{
#' sQTL_DT <- xQTLdownload_sqtlAllAsso(genes=c("MMP7","TP53"), tissue="Lung")
#' }
xQTLdownload_sqtlAllAsso <- function(genes="", geneType="auto", tissue="", clu_names="", clu_geneid_DF=NULL){
  tissueSiteDetail <- clu_name<- tissueSiteDetailId <- gencodeId <- gencodeId_unv <- NULL
  # match tissue:
  if(tissue==""){ stop("tissue can not be null...") }
  tissueDT <- tissueSiteDetailGTExv8[tissueSiteDetail== tissue | tissueSiteDetailId==tissue][1,]
  if(nrow(tissueDT)==0){ stop("tissue not found...")}
  # get clu-gene info:
  if(!is.null(clu_geneid_DF)){
    setDT(clu_geneid_DF)
    clu_geneid<-clu_geneid_DF
  }else{
    clu_geneid <- fread(paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/sQTL_clu_gene2/sQTL_intron_gene_id_",tissueDT$tissueSiteDetailId,".txt"), header=FALSE)
    if(!exists("clu_geneid") || nrow(clu_geneid)==0){ stop("tissue not found...") }
  }
  names(clu_geneid) <- c("clu_name","gencodeId")
  clu_geneid$gencodeId_unv <- unlist(lapply(clu_geneid$gencodeId, function(x){ str_split(x, fixed("."))[[1]][1] }))

  # merge with genes:
  if(geneType=="gencodeId"){
    clu_geneid_sub <- clu_geneid[gencodeId %in% genes |  gencodeId_unv %in% genes]
  }else{
    message("==> Querying genes...")
    genes_info <- xQTLquery_gene(genes)
    clu_geneid_sub <- clu_geneid[gencodeId %in% genes_info$gencodeId]
  }
  # merge with clu_name:
  if(!clu_names==""){
    clu_geneid_sub <- clu_geneid_sub[clu_name %in% clu_names]
  }

  # construct clu infor data.table:
  if(nrow(clu_geneid_sub) == 0){ stop("None gene found...")}
  clu_geneid_sub$url1 <- paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/sQTL_b38_rsid_gene/", tissueDT$tissueSiteDetailId, "/",clu_geneid_sub$clu_name)
  clu_geneid_sub$tmpFilePath <- paste(unlist(lapply(1:nrow(clu_geneid_sub), function(x){tempfile(pattern = "sqtl_")})),"sQTL_",1:nrow(clu_geneid_sub),".txt", sep="")

  # download sQTL by clu:
  message("== Start downloading sQTLs...")
  for(i in 1:nrow(clu_geneid_sub)){
    cat("== Downloading ", clu_geneid_sub[i,]$gencodeId, "-", clu_geneid_sub[i,]$clu_name,"...")
    df <- try(suppressWarnings(utils::download.file(url = clu_geneid_sub[i,]$url1,
                                                    destfile=clu_geneid_sub[i,]$tmpFilePath,
                                                    quiet = TRUE )), silent=TRUE)
    if(!file.exists(clu_geneid_sub[i,]$tmpFilePath)){ cat("    > Failed...") }
    cat("   > Success!")
    message("")
  }

  message("== combine results...")
  sQTL_summary <- data.table()
  for(i in 1:nrow(clu_geneid_sub)){
    if(file.exists(clu_geneid_sub[i,]$tmpFilePath)){
      sQTL_summary_i <- fread(clu_geneid_sub[i,]$tmpFilePath, header=TRUE)
      if(nrow(sQTL_summary_i)>0){
        # cat(i,"-")
        sQTL_summary_i$gencodeId <- clu_geneid_sub[i,]$gencodeId
        sQTL_summary <- rbind(sQTL_summary, sQTL_summary_i)
        file.remove(clu_geneid_sub[i,]$tmpFilePath)
      }
    }
  }
  names(sQTL_summary) <- c("rsid", "clu_name", "pValue", "beta", "se", "gencodeId")
  return(sQTL_summary)
}

#' @title Download summary statistics of xQTL for a specified gene, default:3'aQTL
#'
#' @param genes (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "geneSymbol".
#' @param tissue (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8`
#' @param mRNA_refseq (character) If provided, only the 3'aQTL of mRNA will be downloaded
#' @param mRNA_gene_DF (data.frame) If provided, mRNA-gencode mapping relationship will be loaded from this data.frame.
#' @param type 3'aQTL(default)
#'
#' @return A data.table object of xQTL dataset.
#' @export
#'
#' @examples
#' \donttest{
#' aQTL_DT <- xQTLdownload_xqtlAllAsso(genes=c("MMP7", "EPS15"), tissue="Lung")
#' }
xQTLdownload_xqtlAllAsso <- function(genes="", geneType="geneSymbol",  tissue="", mRNA_refseq="",  mRNA_gene_DF=NULL, type="3'aQTL"){
  tissueSiteDetail <- tissueSiteDetailId <- gencodeId_unv <- geneSymbol <- mRNA<- NULL
  if(type=="3'aQTL"){
    # match tissue:
    tissueDT <- tissueSiteDetailGTExv8[tissueSiteDetail== tissue | tissueSiteDetailId==tissue][1,]
    if(nrow(tissueDT)==0){ stop("tissue not found...")}
    # get clu-gene info:
    if(!is.null(mRNA_gene_DF)){
      setDT(mRNA_gene_DF)
      mRNA_gene<-mRNA_gene_DF
    }else{
      # mRNA_gene <- fread(paste0("https://github.com/dingruofan/exampleData/raw/master/hg38_refseq_gencodeId.txt"), header=TRUE)
      mRNA_gene <- fread("http://bioinfo.szbl.ac.cn/aQTL/for_xQTLbiolinks/hg38_3UTR_anno.txt", header=FALSE)
      if(!exists("mRNA_gene") || nrow(mRNA_gene)==0){ stop("mRNA info download fail...") }
    }
    names(mRNA_gene) <- c("mRNA", "geneSymbol")

    # merge with genes:
    if(geneType=="geneSymbol"){
      mRNA_gene_sub <- mRNA_gene[str_to_lower(geneSymbol) %in% str_to_lower(genes)]
    }else{
      message("==> Querying genes...")
      genes_info <- xQTLquery_gene(genes)
      mRNA_gene_sub <- mRNA_gene[str_to_lower(geneSymbol) %in% str_to_lower(genes_info$geneSymbol)]
    }
    # merge with clu_name:
    if(!mRNA_refseq==""){
      mRNA_gene_sub <- mRNA_gene_sub[str_to_lower(mRNA) %in% str_to_lower(mRNA_refseq)]
    }

    # construct mRNA infor data.table:
    if(nrow(mRNA_gene_sub) == 0){ stop("None gene found...")}
    mRNA_gene_sub$url1 <- paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/aQTL_b38_rsid_gene/", tissueDT$tissueSiteDetailId, "/",mRNA_gene_sub$mRNA)
    mRNA_gene_sub$tmpFilePath <- paste(unlist(lapply(1:nrow(mRNA_gene_sub), function(x){tempfile(pattern = "aqtl_")})),"sQTL_",1:nrow(mRNA_gene_sub),".txt", sep="")

    # download sQTL by clu:
    message("== Start downloading 3'QTLs...")
    for(i in 1:nrow(mRNA_gene_sub)){
      cat("== Downloading ", i,"-", mRNA_gene_sub[i,]$gencodeId, "-", mRNA_gene_sub[i,]$clu_name,"...")
      df <- try(suppressWarnings(utils::download.file(url = mRNA_gene_sub[i,]$url1,
                                                      destfile=mRNA_gene_sub[i,]$tmpFilePath,
                                                      quiet = TRUE )), silent=TRUE)
      if(!file.exists(mRNA_gene_sub[i,]$tmpFilePath)){
        # cat("    > Failed...")
        }else{
        cat("   > Success!")
      }
      message("")
    }

    message("== combine results...")
    aQTL_summary <- data.table()
    for(i in 1:nrow(mRNA_gene_sub)){
      if(file.exists(mRNA_gene_sub[i,]$tmpFilePath)){
        aQTL_summary_i <- fread(mRNA_gene_sub[i,]$tmpFilePath, header=TRUE)
        if(nrow(aQTL_summary_i)>0){
          # cat(i,"-")
          aQTL_summary_i$geneSymbol <- mRNA_gene_sub[i,]$geneSymbol
          aQTL_summary <- rbind(aQTL_summary, aQTL_summary_i)
          file.remove(mRNA_gene_sub[i,]$tmpFilePath)
        }
      }
    }
    names(aQTL_summary) <- c("rsid", "maf","pValue",  "beta", "se", "mRNA", "geneSymbol")
    return(aQTL_summary)
  }
}


#' @title Download all sc-eQTL associations for a specified gene
#'
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "geneSymbol".
#' @param cell_type (character)cell types supported in the list of study_info from 'xQTLquery_scInfo'
#' @param cell_state (character)cell states supported in the list of study_info from 'xQTLquery_scInfo'
#' @param qtl_type (character)QTL types supported in the list of study_info from 'xQTLquery_scInfo'
#' @param study_name (character)study name supported in the list of study_info from 'xQTLquery_scInfo'
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#'  sceQTL_dt <- xQTLdownload_sc(gene="TP53", cell_type = "B Cell", cell_state="-",
#'                  qtl_type="Cell-type-specific eQTL", study_name = "Resztak2022biorxiv")
#' }
xQTLdownload_sc <- function(gene="BIN3",geneType="geneSymbol", cell_type="Astrocytes", cell_state="", qtl_type="Cell-type eQTL", study_name="Bryois2022NN"){
  tmp <- tissueSiteDetail <- genes <- NULL
  study_info <- xQTLquery_scInfo()

  if(geneType=="gencodeId"){
    gene_info <- xQTLquery_gene(gene)
    if(nrow(gene_info)==0){message("invalid gene!");return(data.table())}
    genename  <- gene_info$geneSymbol
  }else if(geneType == "geneSymbol"){
    genename <- gene
  }
  if( !(toupper(cell_type) %in% toupper(unique(study_info$cell_type_name))) ){
    stop(" == cell type |",cell_type,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if( !(toupper(cell_state) %in% toupper(unique(study_info$cell_state))) ){
    stop(" == cell state|",cell_state,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if( !(toupper(qtl_type) %in% toupper(unique(study_info$QTL.type))) ){
    stop(" == QTL type|",qtl_type,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if( !(toupper(study_name) %in% toupper(unique(study_info$study))) ){
    stop(" == study name|",study_name,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if(cell_state==""){cell_state <- "-"}

  url1 <- paste0("http://bioinfo.szbl.ac.cn/scQTLbase_backend/query_xQTLbiolinks?geneName=",genename,"&cellType=",cell_type,"&cellState=",cell_state, "&study=",study_name,"&qtlType=",qtl_type )
  message(url1)
  url1 <- utils::URLencode(url1)
  url1 <- stringr::str_replace_all(url1, fixed("+"), fixed("%2B"))
  url1 <- stringr::str_replace_all(url1, fixed("("), fixed("%28"))
  url1 <- stringr::str_replace_all(url1, fixed(")"), fixed("%29"))
  url1 <- stringr::str_replace_all(url1, fixed(";"), fixed("%3B"))
  print(url1)
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto", retryTimes=11)
  qtl_summary <- data.table::as.data.table(url1GetText2Json$singleTissueEqtl)
  if(nrow(qtl_summary)==0){
    message("No QTL associations were found in ",tissueSiteDetail, " of thess ", length(genes), " genes!")
    message("== Done.")
    return(data.table::data.table())
  }
  return(qtl_summary)
}


#' @title Download significant sc-eQTL associations for a specified gene
#'
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "geneSymbol".
#' @param cell_type (character)cell types supported in the list of study_info from 'xQTLquery_scInfo'
#' @param cell_state (character)cell states supported in the list of study_info from 'xQTLquery_scInfo'
#' @param qtl_type (character)QTL types supported in the list of study_info from 'xQTLquery_scInfo'
#' @param study_name (character)study name supported in the list of study_info from 'xQTLquery_scInfo'
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \donttest{
#'  sceQTL_dt <- xQTLdownload_scSig(gene="CNTNAP2", cell_type = "Definitive Endoderm (defendo)",
#'    cell_state="Response With Pseudotime", qtl_type="Dynamic eQTL", study_name = "cuomo2020NC")
#' }
xQTLdownload_scSig <- function(gene="BIN3",geneType="geneSymbol", cell_type="Astrocytes", cell_state="", qtl_type="Cell-type eQTL", study_name="Bryois2022NN"){
  tmp <- tissueSiteDetail <- genes <- NULL
  study_info <- xQTLquery_scInfo()

  if(geneType=="gencodeId"){
    gene_info <- xQTLquery_gene(gene)
    if(nrow(gene_info)==0){message("invalid gene!");return(data.table())}
    genename  <- gene_info$geneSymbol
  }else if(geneType == "geneSymbol"){
    genename <- gene
  }
  if( !(toupper(cell_type) %in% toupper(unique(study_info$cell_type_name))) ){
    stop(" == cell type |",cell_type,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if( !(toupper(cell_state) %in% toupper(unique(study_info$cell_state))) ){
    stop(" == cell state|",cell_state,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if( !(toupper(qtl_type) %in% toupper(unique(study_info$QTL.type))) ){
    stop(" == QTL type|",qtl_type,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if( !(toupper(study_name) %in% toupper(unique(study_info$study))) ){
    stop(" == study name|",study_name,"| is not in the list. please check using 'xQTLquery_scInfo'")
  }
  if(cell_state==""){cell_state <- "-"}

  url1 <- paste0("http://bioinfo.szbl.ac.cn/scQTLbase_backend/query_xQTLbiolinks_sig?geneName=",genename,"&cellType=",cell_type,"&cellState=",cell_state, "&study=",study_name,"&qtlType=",qtl_type )
  message(url1)
  url1 <- utils::URLencode(url1)
  url1 <- stringr::str_replace_all(url1, fixed("+"), fixed("%2B"))
  url1 <- stringr::str_replace_all(url1, fixed("("), fixed("%28"))
  url1 <- stringr::str_replace_all(url1, fixed(")"), fixed("%29"))
  url1 <- stringr::str_replace_all(url1, fixed(";"), fixed("%3B"))
  # print(url1)
  url1GetText2Json <- fetchContent(url1, method = "download", downloadMethod = "auto", retryTimes=11)
  qtl_summary <- data.table::as.data.table(url1GetText2Json$singleTissueEqtl)
  if(nrow(qtl_summary)==0){
    message("No expression profiles were found in ",tissueSiteDetail, " of thess ", length(genes), " genes!")
    message("== Done.")
    return(data.table::data.table())
  }
  return(qtl_summary)
}

#' @title Export expression object to a specified format
#'
#' @param exp_object expression object derived from `xQTLdownload_exp`
#' @param out_format "to_clusterP", "to_wgcna" and to "to_deseq"
#'
#' @importFrom SummarizedExperiment assay
#'
#' @return A data.frame/data.table object
#' @export
#'
#' @examples
#' \donttest{
#' expProfiles <- xQTLdownload_exp(c("tp53","naDK","SDF4"),
#'                                tissueSiteDetail="Artery - Coronary",
#'                                pathologyNotesCategories=TRUE, toSummarizedExperiment = FALSE)
#' xQTL_export(expProfiles)
#' }
xQTL_export <- function(exp_object, out_format="to_clusterP"){
  if(out_format == "to_clusterP"){
    if(as.character(class(exp_object)) == "RangedSummarizedExperiment"){
      exp_object <- SummarizedExperiment::assay(exp_object)
    }
    gene_info <- xQTLquery_gene(exp_object$gencodeId)[,c("geneSymbol", "gencodeId", "entrezGeneId")]
    return(gene_info)
  }else if(out_format == "to_wgcna"){
    if(as.character(class(exp_object)) == "RangedSummarizedExperiment"){
      exp_object <- SummarizedExperiment::assay(exp_object)
    }
    gene_exp <- as.data.frame(t(exp_object[,which(stringr::str_detect(names(exp_object), "GTEX-"))]))
    names(gene_exp) <- exp_object$geneSymbol
    return(gene_exp)
  }else if(out_format == "to_deseq2"){
    if(as.character(class(exp_object)) == "RangedSummarizedExperiment"){
      exp_object <- SummarizedExperiment::assay(exp_object)
    }
    gene_exp <- as.data.frame(exp_object[,which(stringr::str_detect(names(exp_object), "GTEX-"))])
    names(gene_exp) <- exp_object$geneSymbol
    return(gene_exp)
  }
}



#' @title Download metadata for H3K4me1 and H3K27ac histone QTL (hQTL)
#'
#' @param histone_type (string) One of the histone types: "H3K27AC" or "H3K4ME1".
#' @param cell_type (string) One of the cell types: "monocyte", "neutrophil" or "T cell".
#'
#' @return a data.table object includng all CpG ID
#' @export
#'
#' @examples
#' hqtlmeta <- xQTLdownload_hqtlmeta(histone_type="H3K4ME1", cell_type="T cell")
xQTLdownload_hqtlmeta <- function(histone_type="H3K27AC", cell_type="monocyte"){
  cell_dt <- data.table(cell_types = c("monocyte", "neutrophil", "T cell"), cell_abbr= c("mono", "neut", "tcel"))
  histone_dt <- data.table( histone_types = c("H3K27AC", "H3K4ME1"), histone_abbr=c("K27AC", "K4ME1"))
  cell_abbr <- na.omit(cell_dt[cell_type, on="cell_types"])
  if(nrow(cell_abbr)==0){stop("Invalid cell type, please select from: \"monocyte\", \"neutrophil\" or \"T cell\"")}
  histone_abbr <- na.omit(histone_dt[histone_type, on = "histone_types"])
  if(nrow(histone_abbr)==0){stop("Invalid histone type, please select from: \"H3K27AC\" or \"H3K4ME1\"")}
  #
  histone_cell <- paste0(cell_abbr$cell_abbr, "_", histone_abbr$histone_abbr, "_genes")
  url1 <- paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/haQTL/",histone_cell)
  cat("== Downloading metadata of ",histone_type,"-", cell_type, "for hQTL...")
  meta_hqtl <- fread(url1, header = TRUE)
  meta_hqtl$type <- paste0(cell_abbr$cell_abbr, "_", histone_abbr$histone_abbr)
  if(nrow(meta_hqtl)==0){
    stop("Download fail because of unstable network")
  }
  return(meta_hqtl)
}


#' @title Download summary statistics data of H3K4me1 and H3K27ac histone QTL (hQTL) using a specified location
#'
#' @param phenotype_id phenotype_id that formatted with genome location, like: 9-99773935-99776816, can be obtained using `xQTLdownload_hqtlmeta`
#' @param histone_type (string) One of the histone types: "H3K27AC" or "H3K4ME1".
#' @param cell_type (string) One of the cell types: "monocyte", "neutrophil" or "T cell".
#' @param hqtlmeta A data.table object obtained via `xQTLdownload_hqtlmeta`.
#'
#' @return A data.table object
#' @export
#'
#' @examples
#'  hQTL_dt <- xQTLdownload_hqtl(phenotype_id="10:10458128-10465096",
#'                               histone_type="H3K4ME1", cell_type="T cell")
xQTLdownload_hqtl <- function(phenotype_id="9:99773935-99776816", histone_type="H3K27AC", cell_type="monocyte", hqtlmeta=NULL ){
  #
  phenotype_id_ <- phenotype_id
  if(is.null(hqtlmeta)){
    meta_hqtl <- xQTLdownload_hqtlmeta(histone_type, cell_type)
  }
  meta_hqtl <- na.omit(meta_hqtl[phenotype_id==phenotype_id_,])
  if(nrow(meta_hqtl)==0){stop("Invalid phenotype_id, please use `xQTLdownload_hqtlmeta` to obtain all CpG IDs")}

  # re format:
  phenotype_id_ <- stringr::str_split(phenotype_id_, stringr::regex("[:-]"))[[1]]
  phenotype_id_ <- paste0(paste0(phenotype_id_, collapse = "-"),".txt")
  #
  url1 <- paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/haQTL/", meta_hqtl$type, "/", phenotype_id_)
  hqtl <- fread(url1)
  return(hqtl)
}


#' @title Download metadata of DNA methylation QTL (mQTL)
#'
#' @param tissue_name (String)  One of the tissues: BreastMammaryTissue, ColonTransverse, KidneyCortex, Lung, MuscleSkeletal, Ovary, Prostate, Testis and WholeBlood
#' @return A data.table object
#' @export
#'
#' @examples
#' mQTL_meta<- xQTLdownload_mqtlmeta("Testis")
xQTLdownload_mqtlmeta <- function(tissue_name="BreastMammaryTissue"){
  mQTL_meta <- . <- NULL
  tissues <- c("BreastMammaryTissue", "ColonTransverse", "KidneyCortex", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")
  if(!(tissue_name %in% tissues)){ stop(paste0("Invalid tissue name, must be selected from ", paste0(tissues, collapse = ", "))) }
  #
  url1 <- paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/mQTL/", tissue_name, "_cpg_id.txt")
  temp_x <- tempfile()
  download.file(url1, temp_x)
  if(file.exists(temp_x)){
    mqtl_meta <- fread(temp_x, header = FALSE)
    file.remove(temp_x)
    # mqtl_meta <- fread(url1, header=FALSE)
    names(mqtl_meta) <- "cpg_id"
    mqtl_meta$cpg_id <- str_remove(mqtl_meta$cpg_id, stringr::fixed(".txt"))
    return(mqtl_meta)
  }else{
    stop("download failed")
  }

}


#' @title Download summary statistics data of DNA methylation QTL (mQTL) using CpG ID
#'
#' @param cpg_id phenotype_id like: cg00000236, can be obtained using `xQTLdownload_mqtlmeta`
#' @param tissue_name (String)  One of the tissues: BreastMammaryTissue, ColonTransverse, KidneyCortex, Lung, MuscleSkeletal, Ovary, Prostate, Testis and WholeBlood
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' mQTL_dt <- xQTLdownload_mQTL(cpg_id="cg00000221", tissue_name="Prostate")
xQTLdownload_mQTL <- function(cpg_id="cg00000221",tissue_name="WholeBlood", mQTL_meta=NULL){
  tissues <- c("BreastMammaryTissue", "ColonTransverse", "KidneyCortex", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "WholeBlood")
  if(!(tissue_name %in% tissues)){ stop(paste0("Invalid tissue name, must be selected from ", paste0(tissues, collapse = ", "))) }
  if(!is.null(mQTL_meta)){
    cpg_id_ <- cpg_id
    mQTL_meta <- na.omit(mQTL_meta[cpg_id==cpg_id_,])
    if(nrow(mQTL_meta)==0){ stop("Invalid cpg_id in ", tissue_name)}
  }
  url1 <- paste0("http://bioinfo.szbl.ac.cn/xQTL_biolinks/mQTL/", tissue_name, "/",cpg_id,".txt")
  temp_x <- tempfile()
  download.file(url1, temp_x)

  if(file.exists(temp_x)){
    mqtl <- fread(temp_x, header = TRUE)
    file.remove(temp_x)
    return(mqtl)
  }else{
      stop("Download failed!")
    }
}




