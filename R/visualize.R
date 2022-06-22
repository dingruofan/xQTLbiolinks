#' @title Boxplot of normalized expression among genotypes for eQTL.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @import curl
#' @import jsonlite
#' @return A list containing eQTL detail, expression profile and a ggplot object.
#' @export
#'
#' @examples
#' # EQTL associatons of TP53 in Esophagus - Mucosa:
#' expEqtl <- xQTLvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Liver")
#' expEqtl <- xQTLvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Lung")
#'
#' # EQTL associatons of IRF5:
#' expEqtl<-xQTLvisual_eqtlExp(variantName="rs3778754",gene ="IRF5",tissueSiteDetail="Whole Blood")
xQTLvisual_eqtlExp <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail="", datasetId="gtex_v8" ){
  genoLabels <- normExp <- labelNum <- p <- NULL

  # library(crayon)
  # cat(green(
  #   'I am a green line ' %+%
  #     blue$underline$bold('with a blue substring') %+%
  #     yellow$italic(' that becomes yellow and italicised!\n')
  # ))

  # gene="ENSG00000248746.5"
  # geneType="gencodeId"
  # variantName="chr11_66561248_T_C_b38"
  # variantType="variantId"
  # tissueSiteDetail = "Esophagus - Mucosa"
  # datasetId="gtex_v8"


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
  }else if(gene ==""){
    stop("Parameter \"gene\" can not be null. ")
  }
  # check: variantName
  if(is.null(variantName) ||  any(is.na(variantName)) ){
    stop("Parameter \"variantName\" can not be NULL or NA!")
  }else if(length(variantName)!=1 ||variantName==""){
    stop("Parameter \"variantName\" can not be NULL or NA!")
  }

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
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

  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) || any(is.na(tissueSiteDetail)) || tissueSiteDetail=="" ){
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

  message("== Querying significant eQTL associations from API server:")
  suppressMessages(eqtlInfo <- xQTLdownload_eqtl(gene = gene, variantName = variantName, geneType = geneType, variantType = variantType, tissueSiteDetail = tissueSiteDetail))
  if( !exists("eqtlInfo") || is.null(eqtlInfo) || nrow(eqtlInfo)==0 ){
    stop("No eqtl associations were found for gene [", gene, "] and variant [", variantName,"] in ", tissueSiteDetail, " in ", datasetId,".")
  }else{
    message("   A record of eQTL association was found in ", datasetId, " of gene [", gene,"] and variant [", variantName," in [", tissueSiteDetail,"].")
    message("== Done")
  }

  message("== Querying expression from API server:")
  suppressMessages(eqtlExp <- xQTLdownload_eqtlExp(variantName = eqtlInfo$snpId, gene = eqtlInfo$geneSymbol, tissueSiteDetail = tissueSiteDetail))
  if( !exists("eqtlExp") || is.null(eqtlExp) || nrow(eqtlExp)==0 ){
    stop("No expression profiles were found for gene [", gene, "] in ", tissueSiteDetail, " in ", datasetId,".")
  }else{
    message("   Normalized expression of [",nrow(eqtlExp), "] samples for gene [",gene,"] were obtaioned in ", tissueSiteDetail, " in ", datasetId,".")
    message("== Done")
  }

  eqtlInfo <- cbind(eqtlInfo, data.table::rbindlist(lapply(eqtlInfo$variantId, function(x){var_tmp = stringr::str_split(x,stringr::fixed("_"))[[1]]; data.table::data.table(chrom=var_tmp[1], pos=var_tmp[2], ref=var_tmp[3], alt=var_tmp[4])  })))
  # replace genotypes with ref and alt:
  genoLable <- data.table::data.table(genotypes=0:2, genoLabels = c( ifelse(nchar(eqtlInfo$ref)==1, paste0(rep(eqtlInfo$ref,2),collapse = ""), "Ref"),
                                            ifelse(nchar(eqtlInfo$ref)==1, paste0(c(eqtlInfo$ref, eqtlInfo$alt),collapse = ""), "Het"),
                                            ifelse(nchar(eqtlInfo$ref)==1, paste0(rep(eqtlInfo$alt,2),collapse = ""), "Hom")) )
  genoLable$genoLabels <- factor(genoLable$genoLabels, levels = genoLable$genoLabels)
  genoLable <- merge(eqtlExp, genoLable, by ="genotypes")

  # x axis label:
  genoLableX <- data.table::as.data.table(table(genoLable$genoLabels))
  names(genoLableX) <- c("genoLabels", "Num")
  genoLableX$label <- paste0(genoLableX$genoLabels, "(",genoLableX$Num,")")
  genoLableX <- merge(data.table(genoLabels = levels(genoLable$genoLabels)),genoLableX, by="genoLabels", sort=FALSE)

  # for Pie:
  genoLabelPie <- data.table::data.table(table(genoLable$genoLabels))
  names(genoLabelPie) <- c("genoLabels", "labelNum")
  genoLabelPie$legends <- paste0(genoLabelPie$genoLabels, "(",genoLabelPie$labelNum,")")


  if( requireNamespace("ggplot2") ){
    p<- ggplot( genoLable, aes(x=genoLabels, y=normExp)) +
      geom_violin( aes(fill=genoLabels),width=0.88, trim=FALSE, alpha=0.9, scale="width") +
      geom_boxplot(fill="white", width=0.2,  alpha=0.9)+
      scale_fill_brewer(palette="Dark2") +
      scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
      # labs(title = paste0(ifelse(eqtlInfo$snpId==""|| is.na(eqtlInfo$snpId), eqtlInfo$variantId, eqtlInfo$snpId), "- ", eqtlInfo$geneSymbol) )+
      xlab("Genotypes")+
      ylab("Normalized expression")+
      theme_classic()+
      theme(
        axis.text.x=element_text(size=rel(1.3)),
        axis.text.y = element_text(size=rel(1.3)),
        axis.title = element_text(size=rel(1.3)),
        legend.position = "none",
        plot.title = element_text(hjust=0.5)
      )+
      geom_text(aes(x=ifelse(length(unique(genoLable$genoLabels))==3, 2, 1.5), y=max(genoLable$normExp+1.2), label=paste0("P-value: ",signif(eqtlInfo$pValue, 3)) ))
  }
  return(list(eqtl=eqtlInfo, exp=genoLable, p=p))
}

#' @title Boxplot of normalized expression among genotypes for sQTL.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param phenotypeId A character string. Format like: "chr1:497299:498399:clu_54863:ENSG00000239906.1"
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @import curl
#' @import jsonlite
#' @return A list containing variant detail, expression profile and a ggplot object.
#' @export
#'
#' @examples
#' expSqtl <-xQTLvisual_sqtlExp(variantName="chr11_66561248_T_C_b38",
#'           phenotypeId ="chr11:66348070:66353455:clu_8500:ENSG00000255468.6",
#'           tissueSiteDetail="Skin - Sun Exposed (Lower leg)")
xQTLvisual_sqtlExp <- function(variantName="", phenotypeId="", variantType="auto", tissueSiteDetail="", datasetId="gtex_v8" ){
  genoLabels <- normExp <-geneType<- labelNum <- p <- NULL

  # library(crayon)
  # cat(green(
  #   'I am a green line ' %+%
  #     blue$underline$bold('with a blue substring') %+%
  #     yellow$italic(' that becomes yellow and italicised!\n')
  # ))

  gencodeVersion <- "v26"
  if( datasetId == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }else if( datasetId == "gtex_v7" ){
    gencodeVersion <- "v19"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv7)
  }

  # check: variantName
  if(is.null(variantName) ||  any(is.na(variantName)) ){
    stop("Parameter \"variantName\" can not be NULL or NA!")
  }else if(length(variantName)!=1 ||variantName==""){
    stop("Parameter \"variantName\" can not be NULL or NA!")
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

  # 获得突变信息：
  varInfo <- xQTLquery_varId(variantName=variantName, variantType=variantType,datasetId="gtex_v8" )


  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) || any(is.na(tissueSiteDetail)) || tissueSiteDetail=="" ){
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


  message("== Querying expression from API server:")
  suppressMessages(sqtlExp <- xQTLdownload_sqtlExp(variantName = variantName, variantType = variantType, phenotypeId = phenotypeId, tissueSiteDetail = tissueSiteDetail))
  if( !exists("sqtlExp") || is.null(sqtlExp) || nrow(sqtlExp)==0 ){
    stop("No expression profiles were found for phenotypeId [", phenotypeId, "] in ", tissueSiteDetail, " in ", datasetId,".")
  }else{
    message("   Normalized expression of [",nrow(sqtlExp), "] samples were obtaioned in ", tissueSiteDetail, " in ", datasetId,".")
    message("== Done")
  }

  varInfo <- cbind(varInfo, data.table::rbindlist(lapply(varInfo$variantId, function(x){var_tmp = stringr::str_split(x,stringr::fixed("_"))[[1]]; data.table::data.table(chrom=var_tmp[1], pos=var_tmp[2], ref=var_tmp[3], alt=var_tmp[4])  })))
  # replace genotypes with ref and alt:
  genoLable <- data.table::data.table(genotypes=0:2, genoLabels = c( ifelse(nchar(varInfo$ref)==1, paste0(rep(varInfo$ref,2),collapse = ""), "Ref"),
                                                                     ifelse(nchar(varInfo$ref)==1, paste0(c(varInfo$ref, varInfo$alt),collapse = ""), "Het"),
                                                                     ifelse(nchar(varInfo$ref)==1, paste0(rep(varInfo$alt,2),collapse = ""), "Hom")) )
  genoLable$genoLabels <- factor(genoLable$genoLabels, levels = genoLable$genoLabels)
  genoLable <- merge(sqtlExp, genoLable, by ="genotypes")

  # x axis label:
  genoLableX <- data.table::as.data.table(table(genoLable$genoLabels))
  names(genoLableX) <- c("genoLabels", "Num")
  genoLableX$label <- paste0(genoLableX$genoLabels, "(",genoLableX$Num,")")
  genoLableX <- merge(data.table(genoLabels = levels(genoLable$genoLabels)),genoLableX, by="genoLabels", sort=FALSE)

  # for Pie:
  genoLabelPie <- data.table::data.table(table(genoLable$genoLabels))
  names(genoLabelPie) <- c("genoLabels", "labelNum")
  genoLabelPie$legends <- paste0(genoLabelPie$genoLabels, "(",genoLabelPie$labelNum,")")

  # retrieve sQTL detail:
  P_geneName <- stringr::str_split(phenotypeId, ":")[[1]][5]
  try(suppressMessages(sqtlInfo <- xQTLdownload_sqtlSig(variantName = varInfo$variantId, genes = P_geneName, tissueSiteDetail=tissueSiteDetail)), silent = TRUE)
  sqtlInfo <- sqtlInfo[phenotypeId==phenotypeId]
  if(exists("sqtlInfo") & !is.null(sqtlInfo) & nrow(sqtlInfo)==1 ){
    message("Significant sQTL association found!")
    labelPvalue <- paste0("P-value: ",signif(sqtlInfo$pValue, 3))
  }else{
    labelPvalue <- ""
  }


  if( requireNamespace("ggplot2") ){
    p<- ggplot( genoLable, aes(x=genoLabels, y=normExp)) +
      geom_violin( aes(fill=genoLabels),width=0.88, trim=FALSE, alpha=0.9, scale="width") +
      geom_boxplot(fill="white", width=0.2,  alpha=0.9)+
      scale_fill_brewer(palette="Dark2") +
      scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
      # labs(title = paste0(ifelse(varInfo$snpId==""|| is.na(varInfo$snpId), varInfo$variantId, varInfo$snpId), "- ", phenotypeId, " (",tissueSiteDetail,")") )+
      xlab("Genotypes")+
      ylab("Normalized expression")+
      theme_classic()+
      theme(
        axis.text=element_text(size=rel(1.3)),
        axis.title = element_text(size=rel(1.3)),
        legend.position = "none",
        plot.title = element_text(hjust=0.5)
      )+
      geom_text(aes(x=ifelse(length(unique(genoLable$genoLabels))==3, 2, 1.5), y=max(genoLable$normExp+1.2), label= labelPvalue))
  }

  return(list(varInfo=varInfo, exp=genoLable, p=p))
}



#' @title Locuszoom plot for visualizing regional signals relative to genomic position with a file of summary statistics
#' @description
#' This function is rebuilt from `locuscompare.R` (https://github.com/boxiangliu/locuscomparer/blob/master/R/locuscompare.R).
#' @param DF A data.frame or a data.table object. Four columns are required (arbitrary column names is supported):
#'
#'  `Col 1`. "snps" (character), , using an rsID (e.g. "rs11966562");
#'
#'  `Col 2`. "chromosome" (character), one of the chromosome from chr1-chr22;
#'
#'  `Col 3`. "postion" (integer), genome position of snp.
#'
#'  `Col 4`. "P-value" (numeric).
#' @param highlightSnp Default is the snp that with lowest p-value.
#' @param population One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @param posRange Genome range that you want to visualize (e.g. "chr6:3e7-7e7"). Default is the region that covers all snps.
#' @param legend (boolean, optional) Should the legend be shown? Default: TRUE.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param snpLD A data.frame of LD matirx. Default is null.
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @importFrom  cowplot ggdraw draw_label
#' @import utils
#' @return A list containing data.table and ggplot object.
#' @export
#'
#' @examples
#' # For GWAS dataset:
#' library(data.table)
#' # load data:
#' gwasDF <- fread("https://gitee.com/stronghoney/exampleData/raw/master/gwasChr6Sub4.txt")
#' xQTLvisual_locusZoom(gwasDF)
#' # Zoom in:
#' xQTLvisual_locusZoom(gwasDF, posRange="chr6:4.7e7-4.8e7", population ="EUR")
#'
#' # For eQTL of a gene of interest:
#' eqtlAsso <- xQTLdownload_eqtlAllAsso("RP11-385F7.1", tissueLabel = "Brain - Cortex",
#'                                      withB37VariantId=FALSE)
#' xQTLvisual_locusZoom(eqtlAsso[,c("snpId", "chrom", "pos", "pValue")], highlightSnp="rs4711878" )
#' # Zoom in:
#' xQTLvisual_locusZoom(eqtlAsso[,c("snpId", "chrom", "pos", "pValue")], highlightSnp="rs4711878",
#'                      posRange="chr6:47.3e6-47.9e6")
xQTLvisual_locusZoom <- function( DF , highlightSnp="", population="EUR", posRange="", legend = TRUE, legend_position = c('topright','bottomright','topleft'),  snpLD=NULL){
  snpId <- pos <- pValue <- logP <- pointShape<- NULL
  chrom <- x <- y<- RS_Number <- R2 <- SNP_B <- r2Cut <-genome<- .<-NULL
  # highlightSnp=""
  # population="EUR"
  # posRange="chr6:46488310-48376712"
  # token="9246d2db7917"
  # windowSize=1e6

  DF <- na.omit(DF[,1:4])
  names(DF) <- c("snpId", "chrom", "pos", "pValue")
  data.table::setDT(DF)
  DF$pos <- as.integer(DF$pos)
  # DF[,chrom:=.(ifelse(stringr::str_detect(chrom,"^chr"), chrom,paste("chr",chrom,sep="")))]
  DF$chrom <- ifelse(stringr::str_detect(DF$chrom,"^chr"), DF$chrom, paste("chr",DF$chrom,sep=""))
  # chrome check:
  P_chrom <- unique(DF$chrom)
  if( length(P_chrom) !=1 || !(P_chrom %in% paste0("chr",c(1:22,"x"))) ){
    stop("SNPs must located in the same chromosome! [", paste(P_chrom, collapse = ","),"] are detected!")
  }

  setDT(DF)
  # remove snp without dbSNP ID:
  DF <- DF[stringr::str_detect(snpId, stringr::regex("^rs")),]
  # retain snps in range:
  if(posRange!=""){
    posRangeSplit <- stringr::str_split(posRange, stringr::regex(":|-"))[[1]]
    DF <-DF[pos>min(as.integer(posRangeSplit[2:3])) & pos< max(as.integer(posRangeSplit[2:3]))]
  }
  if(nrow(DF)<2){
    stop("Please enlarge the genome range!")
  }

  hSnpCount <- 1
  # highligth SNP:
  if(highlightSnp ==""){
    highlightSnp <- DF[order(pValue)][hSnpCount,]$snpId
  }
  message("== Highlighted SNP: [",highlightSnp,"]...")

  # LD info:
  if( is.null(snpLD) ){
    message("== Retrieve LD information of SNP: [",highlightSnp,"]...")
    try(snpLD <- retrieveLD(DF[order(pValue)][hSnpCount,]$chrom, highlightSnp, population))
    # try(snpLD <- retrieveLD_LDproxy(highlightSnp,population = population, windowSize = windowSize, genomeVersion = genomeVersion, token = token) )
    data.table::setDT(snpLD)
  }
  # snpLD <- snpLD[,.(SNP_A=highlightSnp, SNP_B=RS_Number, R2)]

  # Set LD SNP color:
  if( nrow(snpLD)>0 ){
    # 由于多个人种存在重复LD，所以取均值：
    snpLD <- snpLD[,.(R2=mean(R2)),by=c("SNP_A","SNP_B")]
    # snpLD$colorP = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('#636363','#7fcdbb','darkgreen','#feb24c','gold'), include.lowest=TRUE))
    snpLD$r2Cut = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE))
    # snpLD$sizeP = as.character(cut(snpLD$R2,breaks=c(0,0.8, 0.9,1), labels=c(1,1.01,1.1), include.lowest=TRUE))
  }else{
    message("No LD information of [",highlightSnp,"] was detected.")
    snpLD <- data.table(SNP_A=character(0), SNP_B =character(0),R2=numeric(0), color=character(0),r2Cut=character(0) )
  }

  # set color:
  DF$logP <- (-log(DF$pValue, 10))
  DF <- merge(DF, snpLD[,.(snpId=SNP_B, r2Cut)], by ="snpId", all.x=TRUE, sort=FALSE)
  DF[is.na(r2Cut),"r2Cut"]<- "(0.0-0.2]"
  DF[snpId==highlightSnp,"r2Cut"]<- "(0.8-1.0]"

  # set size:
  DF$sizeP <- "small"
  DF[logP<max(DF$logP), "sizeP"] <- "small"
  DF[logP>=2 & logP< max(DF$logP), "sizeP"] <- "middle"
  DF[logP == max(DF$logP), "sizeP"] <- "large"
  DF[snpId == highlightSnp, "sizeP"] <- "most"

  # set color:
  colorDT <- data.table( r2Cut = as.character(cut(c(0.2,0.4,0.6,0.8,1),breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE)),
                         pointFill= c('blue4','skyblue','darkgreen','orange','red'),
                         # pointFill = c("#9C8B88", "#e09351", "#df7e66", "#b75347", "#096CFD"),
                         pointColor = c('black','black','black','black','black'),
                         pointSize = c(1,1,2,2,2.5))
  colorDT <- merge(colorDT, unique(DF[,.(r2Cut)]), by="r2Cut",all.x=TRUE)[order(-r2Cut)]
  # set shape:
  DF$pointShape <- "normal"
  DF[snpId==highlightSnp,"pointShape"] <- "highlight"
  DF$pointShape <- as.factor(DF$pointShape)
  DF <- DF[order(r2Cut, logP)]
  DF <- rbind( DF[snpId!=highlightSnp ], DF[snpId==highlightSnp, ] )

  # title:
  plotTitle <- paste0( ifelse(stringr::str_detect(P_chrom, stringr::regex("^chr")),P_chrom, paste0("chr", P_chrom)),  ":", paste0(range(DF$pos), collapse = "-") )

  # ylab and unit:
  posUnit <- "Bb"
  if( any(range(DF$pos)>10^6)){
    DF$pos <- DF$pos/10^6
    posUnit <- "Mb"
  }else if( all(range(DF$pos)<10^6) && all(range(DF$pos)>10^3) ){
    DF$pos <- DF$pos/10^3
    posUnit <- "Kb"
  }else{
    posUnit <- "Bb"
  }
  yLab <- expression(-log["10"]("Pvalue"))

  # xlab:
  xLab <- paste0(ifelse(stringr::str_detect(P_chrom, stringr::regex("^chr")),P_chrom, paste0("chr", P_chrom))," (",posUnit,")")

  p <- ggplot(DF)+
    geom_point(aes(x=pos, y=logP, fill=r2Cut,  size=pointShape, shape=pointShape), color="black")+
    scale_size_manual(breaks = c('normal', "highlight"), values =  c(3,3.5)  )+
    scale_shape_manual(breaks = c('normal', "highlight"), values =  c(21,23) )+
    scale_y_continuous(breaks = seq(0,floor(max(DF$logP))), labels = seq(0,floor(max(DF$logP))), limits = c(0, max(DF$logP)+0.4))+
    # scale_color_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointColor) +
    scale_fill_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointFill) +
    # geom_text(aes(x=pos, y=logP, label=snpId ))+
    ggrepel::geom_label_repel(data=DF[snpId==highlightSnp,], aes(x=pos, y=logP, label=snpId) )+
    # labs(title = plotTitle )+
    xlab( xLab )+
    ylab( yLab )+
    theme_classic()+
    theme(axis.text=element_text(size=rel(1.5), color = "black"),
          axis.title=element_text(size=rel(1.4), color = "black"),
          plot.title = element_text(hjust=0.5),
          legend.position = "none"
    )+ guides( shape="none", color="none", size="none", fill = "none" )

  if (legend == TRUE) {
    legend_position = match.arg(legend_position)
    if (legend_position == 'bottomright'){
      legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
    } else if (legend_position == 'topright'){
      legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
    } else {
      legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
    }

    p <- cowplot::ggdraw(p) +
      geom_rect(data = legend_box,
                aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                color = "black",
                fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
      draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
      draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
      draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
      draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
      draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
  }

  # print(p)
  return(p)
}

#' @title Dotplot of comparing regional signals between GWAS and xQTL.
#' @description This function is rebuilt from `locuscompare.R` (https://github.com/boxiangliu/locuscomparer/blob/master/R/locuscompare.R).
#' @param eqtlDF A data.frame or data.table with two columns: dbSNP id and p-value.
#' @param gwasDF A data.frame or data.table with two columns: dbSNP id and p-value.
#' @param highlightSnp Default is the snp that is farthest from the origin of the coordinates.
#' @param population One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.#' @param token LDlink provided user token, default = NULL, register for token at https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param legend (boolean, optional) Should the legend be shown? Default: TRUE.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param snpLD A data.frame object of LD matrix. Default is null.
#' @import data.table
#' @import ggplot2
#' @import stringr
#' @import ggrepel
#' @importFrom  cowplot ggdraw draw_label
#' @return A ggplot object.
#' @export
#'
#' @examples
#' library(data.table)
#' # load data:
#' eqtlDF <-fread("https://gitee.com/stronghoney/exampleData/raw/master/eqtl/eqtlAsso1.txt")
#' gwasDF <-fread("https://gitee.com/stronghoney/exampleData/raw/master/gwas/AD/gwasChr6Sub3.txt")
#' # visualize:
#' xQTLvisual_locusCompare( eqtlDF, gwasDF, legend_position="topleft")
xQTLvisual_locusCompare <- function(eqtlDF, gwasDF, highlightSnp="", population="EUR", legend = TRUE, legend_position = c('topright','bottomright','topleft'),  snpLD=NULL ){
  x <- y<- genomeVersion <- NULL

  pValue <- snpId <- distance <- logP.gwas <- logP.eqtl <- NULL
  RS_Number <- R2 <- SNP_B <- r2Cut <- pointShape<- .<-NULL
  eqtlDF <- na.omit(eqtlDF[,1:2])
  gwasDF <- na.omit(gwasDF[,1:2])

  data.table::setDT(eqtlDF)
  data.table::setDT(gwasDF)
  colnames(eqtlDF) <- c("snpId","pValue")
  colnames(gwasDF) <- c("snpId","pValue")
  # remove duplicates:
  eqtlDF <- eqtlDF[order(pValue)][!duplicated(snpId)]
  gwasDF <- gwasDF[order(pValue)][!duplicated(snpId)]

  # remove snp without dbSNP ID:
  eqtlDF <- eqtlDF[stringr::str_detect(snpId, stringr::regex("^rs")),]
  gwasDF <- gwasDF[stringr::str_detect(snpId, stringr::regex("^rs")),]

  # mege:
  DF <- merge(eqtlDF, gwasDF, by="snpId", suffixes = c(".eqtl", ".gwas"))
  if(nrow(DF)<2){
    stop("No shared variants detected, please enlarge the genome range.")
  }

  # log:
  DF$logP.eqtl <- (-log(DF$pValue.eqtl, 10))
  DF$logP.gwas <- (-log(DF$pValue.gwas, 10))
  DF$distance <- sqrt(DF$logP.gwas^2+DF$logP.eqtl^2)

  hSnpCount <- 1
  # highligth SNP:
  if(highlightSnp ==""){
    highlightSnp <- DF[order(-distance)][hSnpCount,]$snpId
    message("== Highlighted SNP: [",highlightSnp,"]...")
  }

  # LD info:
  if( is.null(snpLD) ){
    highlightSnpInfo <- xQTLquery_varId(highlightSnp)
    if(nrow(highlightSnpInfo)==0){
      stop(" Highlighted SNP [", highlightSnp,"] is not detected in GTEx, please set the correct highlightSnp.")
    }
    message(" == Done.")
    message("== Retrieve LD information of SNP: [",highlightSnp,"]...")
    try(snpLD <- retrieveLD(highlightSnpInfo$chromosome, highlightSnp, population))
    # try(snpLD <- retrieveLD_LDproxy(highlightSnp,population = population, windowSize = windowSize, genomeVersion = genomeVersion, token = token) )
    data.table::setDT(snpLD)
  }
  # snpLD <- snpLD[,.(SNP_A=highlightSnp, SNP_B=RS_Number, R2)]

  # Set LD SNP color:
  if( nrow(snpLD)>0 ){
    # 由于多个人种存在重复LD，所以取均值：
    snpLD <- snpLD[,.(R2=mean(R2)),by=c("SNP_A","SNP_B")]
    # snpLD$colorP = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('#636363','#7fcdbb','darkgreen','#feb24c','gold'), include.lowest=TRUE))
    snpLD$r2Cut = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE))
    # snpLD$sizeP = as.character(cut(snpLD$R2,breaks=c(0,0.8, 0.9,1), labels=c(1,1.01,1.1), include.lowest=TRUE))
  }else{
    message("No LD information of SNP [",highlightSnp,"] was detected.")
    snpLD <- data.table::data.table(SNP_A=character(0), SNP_B =character(0),R2=numeric(0), color=character(0),r2Cut=character(0) )
  }
  message("== Done.")

  # set color:
  gwas_eqtl <- merge(DF, snpLD[,.(snpId=SNP_B, r2Cut)], by ="snpId", all.x=TRUE, sort=FALSE)
  gwas_eqtl[is.na(r2Cut),"r2Cut"]<- "(0.0-0.2]"
  gwas_eqtl[snpId==highlightSnp,"r2Cut"] <- "(0.8-1.0]"
  # set color:
  colorDT <- data.table::data.table( r2Cut = as.character(cut(c(0.2,0.4,0.6,0.8,1),breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE)),
                                     pointFill= c('blue4','skyblue','darkgreen','orange','red'),
                                     pointColor = c('black','black','black','black','black')
                                     )
  colorDT <- merge(colorDT, unique(gwas_eqtl[,.(r2Cut)]), by="r2Cut",all.x=TRUE)[order(-r2Cut)]
  # gwas_eqtl <- merge( gwas_eqtl, colorDT, by="r2Cut")
  # set shape:
  gwas_eqtl$pointShape <- "normal"
  gwas_eqtl[snpId==highlightSnp,"pointShape"] <- "highlight"
  gwas_eqtl$pointShape <- as.factor(gwas_eqtl$pointShape)
  gwas_eqtl <- gwas_eqtl[order(r2Cut, logP.gwas, logP.eqtl)]
  gwas_eqtl <- rbind( gwas_eqtl[snpId!=highlightSnp ], gwas_eqtl[snpId==highlightSnp, ] )

  # title:
  plotTitle <- paste0(nrow(gwas_eqtl)," snps")

  # Xlab and ylab:
  xLab <- expression(-log["10"]("p-value") (eQTL))
  yLab <- expression(-log["10"]("p-value") (GWAS))

  if( requireNamespace("ggplot2") ){
    p <- ggplot(gwas_eqtl)+
      geom_point(aes(x=logP.eqtl, y=logP.gwas, fill=r2Cut, color=r2Cut, shape=pointShape, size=pointShape))+
      scale_fill_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointFill)+
      scale_color_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointColor)+
      scale_shape_manual("Highlight",breaks = c('normal', "highlight"), values =  c(21,23) )+
      scale_size_manual("Highlight",breaks = c('normal', "highlight"), values =  c(3,3.5) )+
      # geom_text(aes(x=pos, y=logP, label=snpId ))+
      geom_label_repel(data=gwas_eqtl[snpId==highlightSnp,], aes(x=logP.eqtl, y=logP.gwas, label=snpId) )+
      # labs(title = plotTitle )+
      xlab( xLab )+
      ylab( yLab )+
      theme_classic()+
      theme(axis.text=element_text(size=rel(1.5), color = "black"),
            axis.title=element_text(size=rel(1.5),color = "black"),
            plot.title = element_text(hjust=0.5),
            legend.title = element_text(size=rel(1.3)),
            legend.text = element_text(size=rel(1.2))
      )+ guides( shape="none", color="none", size="none", fill = "none" )

    if (legend == TRUE) {
      legend_position = match.arg(legend_position)
      if (legend_position == 'bottomright'){
        legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
      } else if (legend_position == 'topright'){
        legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
      } else {
        legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
      }

      p = cowplot::ggdraw(p) +
        geom_rect(data = legend_box,
                  aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                  color = "black",
                  fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
        draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
        draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
        draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
        draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
        draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 12)
    }
    # print(p)
  }
  return(p)
}


#
#' @title Generate a combined figure including locuszoom and locuscompare plot object.
#' @description
#' This function is rebuilt from `locuscompare.R` (https://github.com/boxiangliu/locuscomparer/blob/master/R/locuscompare.R).
#' @param gwasEqtldata A data.frame or a data.table that including signals from both GWAS and eQTL. Five columns are required (arbitrary column names is supported):
#'
#'  `Col 1`. "snps" (character), using an rsID (e.g. "rs11966562").
#'
#'  `Col 2`. "chromosome" (character), one of the chromosome from chr1-chr22.
#'
#'  `Col 3`. "postion" (integer), genome position of snp.
#'
#'  `Col 4`. "P-value" (numeric) of GWAS signals.
#'
#'  `Col 5`. "P-value" (numeric) of eQTL signals.
#'
#' @param posRange Genome range that you want to visualize (e.g. "chr6:3e7-7e7"). Default is the region that covers all snps.
#' @param population One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @param highlightSnp Default is the snp that with lowest p-value.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param snpLD A data.frame object of LD matrix. Default is null.
#' @import data.table
#' @import stringr
#' @return A ggplot object.
#' @export
#'
#' @examples
#' # load data:
#' u1 <-"https://raw.githubusercontent.com/dingruofan/exampleData/master/gwas/AD/gwasEqtldata.txt"
#' gwasEqtldata <- data.table::fread(u1)
#' xQTLvisual_locusCombine(gwasEqtldata, highlightSnp="rs13120565")
xQTLvisual_locusCombine <- function(gwasEqtldata, posRange="", population="EUR", highlightSnp="", legend_position="bottomright", snpLD=NULL){
  position <- distance <- rsid <- pValue.eqtl <- pValue.gwas <- chrom <- NULL
  . <- NULL

  gwasEqtldata <- gwasEqtldata[,1:5]
  data.table::setDT(gwasEqtldata)
  names(gwasEqtldata) <- c("rsid", "chrom", "position", "pValue.gwas", "pValue.eqtl")

  # refine chrom:
  # gwasEqtldata[,chrom:=.(ifelse(stringr::str_detect(chrom,"^chr"), chrom, paste("chr",chrom,sep="")))]
  gwasEqtldata$chrom <- ifelse(stringr::str_detect(gwasEqtldata$chrom,"^chr"), gwasEqtldata$chrom, paste("chr",gwasEqtldata$chrom,sep=""))

  # retain snps in range:
  if(posRange!=""){
    posRangeSplit <- stringr::str_split(posRange, stringr::regex(":|-"))[[1]]
    gwasEqtldata <-gwasEqtldata[position>min(as.integer(posRangeSplit[2:3])) & position< max(as.integer(posRangeSplit[2:3]))]
  }
  if(nrow(gwasEqtldata)<2){
    stop("No variant detected in this range, please enlarge the genome range!")
  }

  DF <- data.table::copy(gwasEqtldata)
  DF$logP.eqtl <- (-log(DF$pValue.eqtl, 10))
  DF$logP.gwas <- (-log(DF$pValue.gwas, 10))
  DF$distance <- sqrt(DF$logP.gwas^2+DF$logP.eqtl^2)
  hSnpCount <- 1
  if(highlightSnp ==""){
    highlightSnp <- DF[order(-distance)][hSnpCount,]$rsid
    message("== Highlighted SNP: [",highlightSnp,"]...")
    highlightSnpInfo <- xQTLquery_varId(highlightSnp)
  }else{
    message("== Highlighted SNP: [",highlightSnp,"]")
    highlightSnpInfo <- xQTLquery_varId(highlightSnp)
  }

  # LD info:
  if( is.null(snpLD) ){
    message("== Retrieve LD information of SNP: [",highlightSnp,"]...")
    try( snpLD <- retrieveLD(DF[1,]$chrom, highlightSnp, population) )
    # try(snpLD <- retrieveLD_LDproxy(highlightSnp,population = population, windowSize = windowSize, genomeVersion = genomeVersion, token = token) )
    data.table::setDT(snpLD)
  }

  message("Start plotting locuscomappre...")
  p_scatter<- xQTLvisual_locusCompare(gwasEqtldata[,.(rsid, pValue.eqtl)], gwasEqtldata[,.(rsid, pValue.gwas)],
                                      highlightSnp=highlightSnp, population = population, legend_position = legend_position, snpLD = snpLD)
  message("Start plotting locuszoom for gwas...")
  p_gwas <- xQTLvisual_locusZoom(gwasEqtldata[,.(rsid, chrom, position, pValue.gwas)], legend=FALSE, highlightSnp = highlightSnp, population = population, snpLD = snpLD)
  message("Start plotting locuszoom for eQTL... ")
  p_eqtl <- xQTLvisual_locusZoom(gwasEqtldata[,.(rsid, chrom, position, pValue.eqtl)], legend=FALSE, highlightSnp = highlightSnp, population = population, snpLD = snpLD)

  p_gwas_new = p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  p_gwas_eqtl = cowplot::plot_grid(p_gwas_new, p_eqtl, align = "v", nrow = 2, rel_heights=c(0.8,1))
  p_combined = cowplot::plot_grid(p_scatter, p_gwas_eqtl)
  return(p_combined)
}


#' @title Density plot of expression profiles of the gene
#' @param genes (character string or a character vector) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import ggpubr
#' @importFrom SummarizedExperiment assay colData
#' @return A ggplot object.
#' @export
#'
#' @examples
#' genes <- c("FNDC8", "S100Z", "AQP6", "AMOT", "C3orf38", "FOXL1", "COX11",
#'            "FCN3", "DDX58", "CFI", "MS4A18", "NUDT13", "HOXA4", "VSX1")
#' xQTLvisual_genesExp(genes, tissueSiteDetail="Lung")
#'
#' genes <- c("ENSG00000073598.5","ENSG00000171643.13","ENSG00000086159.12","ENSG00000126016.15",
#'            "ENSG00000179021.9","ENSG00000176678.5","ENSG00000166260.10","ENSG00000142748.12",
#'            "ENSG00000107201.9","ENSG00000205403.12","ENSG00000214782.7","ENSG00000166321.13",
#'            "ENSG00000197576.13","ENSG00000100987.14")
#' xQTLvisual_genesExp(genes, geneType="gencodeId", tissueSiteDetail="Liver")
xQTLvisual_genesExp <- function(genes, geneType="auto", tissueSiteDetail = "", datasetId="gtex_v8"){
  `..density..`<-geneSymbol <- NULL

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  expProfiles <- xQTLdownload_exp(genes=genes, geneType = geneType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId, toSummarizedExperiment=TRUE)
  expData <- as.data.table(cbind( data.table(geneSymbol=rownames(expProfiles)), SummarizedExperiment::assay(expProfiles) ))
  expData1 <- melt(expData, id.vars="geneSymbol", variable.name = "sampleId", value.name="exp")
  p <- ggplot( expData1, aes(x = log(exp+1,10), y = reorder(geneSymbol, -exp, median), fill = ..density..))+
    ggridges::geom_density_ridges_gradient( gradient_lwd = 1, scale = 1.4, rel_min_height = 0.05, size = 0.3) +
    # scale_fill_gradientn( colours = colorRampPalette(c("white", "blue", "red"))(27) )+
    viridis::scale_fill_viridis(name = "WCAE per SOC", option = "C")+
    xlab(expression("Gene expression -log"["10"]("TPM+1")))+
    ylab("Gene symbol")+
    theme_bw()+
    theme( legend.title = element_blank(),
           axis.text.x=element_text(size=rel(1.1),face="bold"),
           axis.text.y = element_text(size=rel(1.1),face="bold"),
           axis.title.x = element_text(size=rel(1.1),face="bold"),
           axis.title.y = element_blank())
  print(p)
  return(p)
}


#' @title Scatter plot for showing the correlation of two genes’ expression.
#' @param gene2 Gene symbol or gencode ID of two genes. Default: gene symbol.
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param groupBy Default:sex, can be choosen from pathologyNotesCategories, like: pathologyNotesCategories.mastopathy, pathologyNotesCategories.mastopathy.metaplasia.
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import stringr
#' @import ggpubr
#' @return A ggplot object.
#' @export
#'
#' @examples
#' gene2 = c("AMOT", "HOXA4")
#' xQTLvisual_geneCorr(gene2,tissueSiteDetail="Liver")
#' xQTLvisual_geneCorr(gene2,groupBy="pathologyNotesCategories.congestion",tissueSiteDetail="Lung")
xQTLvisual_geneCorr <- function(gene2="", geneType="auto", tissueSiteDetail = "", groupBy="sex", datasetId="gtex_v8"){
  geneSymbol <- NULL

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(gene2, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  #
  expProfiles <- xQTLdownload_exp(genes=gene2, geneType = geneType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId, toSummarizedExperiment=TRUE, pathologyNotesCategories = TRUE)
  expData <- as.data.table(cbind( data.table(geneSymbol=rownames(expProfiles)), assay(expProfiles) ))
  expData2 <- as.data.frame(t(expData[geneSymbol %in% gene2][,-c("geneSymbol")]))
  colnames(expData2) <- gene2


  corP <- cor.test(unlist(expData2[,1]), unlist(expData2[,2]))
  message("pearson correlation coefficient: ", corP$estimate," Pvalue: ",corP$p.value)

  sampleInfo <- as.data.table(colData(expProfiles))
  expData2 <- cbind(sampleInfo[rownames(expData2), on="sampleId"][,c("sampleId", groupBy),with=FALSE], expData2)
  expData2 <- na.omit(expData2)

  p <- ggpubr::ggscatterhist(expData2, x=gene2[1], y=gene2[2],
                shape = 21, color = groupBy, fill=groupBy,
                margin.plot="density",
                margin.params = list(fill=groupBy, color="black", size=0.2),
                legend = c(0.9,0.15),
                ggtheme = theme_minimal())
  print(p)
  return(p)
}



#' @title Box plot with jittered points for showing number and significance of eQTL associations
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import PupillometryR
#' @return A ggplot object.
#' @export
#'
#' @examples
#' xQTLvisual_eqtl("KIF15")
#' xQTLvisual_eqtl("MLH1")
xQTLvisual_eqtl <- function(gene, geneType="auto", datasetId = "gtex_v8" ){
  variantId <- tissueSiteDetail <- pValue <- logP <- NULL
  . <- NULL
  # gene="KIF15"
  if( datasetId=="gtex_v8" ){
    gencodeVersion="v26"
  }else{
    gencodeVersion="v19"
  }

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  geneInfo <- xQTLquery_gene(gene, geneType = geneType, gencodeVersion = gencodeVersion )
  geneEqtl <- xQTLdownload_eqtlSig(genes=geneInfo$geneSymbol, datasetId=datasetId)
  geneEqtlSub <- geneEqtl[,.(variantId, tissueSiteDetail, pValue)]
  geneEqtlSub$logP <- -log(geneEqtlSub$pValue, 10)
  setDF(geneEqtlSub)
  p<- ggplot(geneEqtlSub, aes(x=reorder(tissueSiteDetail, -logP, median),y=logP))+
    PupillometryR::geom_flat_violin(data=geneEqtlSub, mapping=aes(fill=tissueSiteDetail), position=position_nudge(x=0.25), color="black", scale = "width")+
    geom_jitter(aes(color=tissueSiteDetail), width = 0.1)+
    geom_boxplot(width=0.15, position = position_nudge(x=0.25), fill="white", size=0.1)+
    coord_flip() +
    ylab(expression(-log["10"]("Pvalue")))+
    xlab("") +
    theme_bw() +
    theme(
      axis.text.x=element_text(size=rel(1.2),face="bold"),
      axis.text.y = element_text(size=rel(1.2),face="bold"),
      axis.title.x = element_text(size=rel(1.3),face="bold"),
      axis.title.y = element_blank()
    )+
    guides(fill="none", color="none")
  print(p)
  return(p)
}


#' @title Violin plot of distribution of the gene expression profiles among multiple tissues.
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissues A character string or a vector. "All" (default) means that all tissues is included.
#' @param datasetId (character) options: "gtex_v8" (default), "gtex_v7".
#' @param log10y Display values of expression in log scale. Default: FALSE.
#' @param toTissueSite TRUE or FALSE, display all subtissues or tissue Site. Default: TURE.
#'
#' @return A list containing expression profile and a ggplot object.
#' @export
#'
#' @examples
#' # Display gene expression in all tissues.
#' # geneExpTissues <- xQTLvisual_geneExpTissues("TP53")
#'
#' # Display gene expression in specified tissues.
#' geneExpTissues <- xQTLvisual_geneExpTissues("TP53", tissues=c("Lung", "Brain","Ovary"))
#'
#' # Display gene expression in log scale in specified tissues.
#' geneExpTissues <- xQTLvisual_geneExpTissues("TP53", tissues="Blood Vessel", log10y=TRUE)
#'
#' # Display gene expression in whole tissue.
#' geneExpTissues <- xQTLvisual_geneExpTissues("TP53", tissues="Blood Vessel", toTissueSite=TRUE)
xQTLvisual_geneExpTissues <- function(gene="", geneType="auto", tissues="All", datasetId="gtex_v8", log10y=FALSE, toTissueSite=FALSE){
  colorHex <- tissueSite <- expTPM <- NULL
  .<-NULL

  if(datasetId == "gtex_v8"){
    tissueSiteDetail <- data.table::copy(tissueSiteDetailGTExv8)
  }else if(datasetId == "gtex_v7"){
    tissueSiteDetail <- data.table::copy(tissueSiteDetailGTExv7)
  }else{
    stop("Please choose the right datasetId!")
  }
  tissueSiteDetail <- tissueSiteDetail[order(tissueSite, tissueSiteDetail)]

  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(gene, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  expProfiles <- data.table()
  if( length(tissues) ==1 && tissues=="All" ){
    tissueSiteDetail_ <- tissueSiteDetail$tissueSiteDetail
  }else{
    tissueSiteDetail_ <- rbind(tissueSiteDetail[tissueSiteDetail %in% tissues], tissueSiteDetail[tissueSite %in% tissues])$tissueSiteDetail
  }
  message("== Start fetching expression profiles of gene [",gene,"] in following tissues...")
  for( tt in 1:length(tissueSiteDetail_)){
    message("   ", tt, " | ",  tissueSiteDetail_[tt])
  }
  message("== This may take A few minutes...", format(Sys.time(), " | %Y-%b-%d %H:%M:%S "))
  for(t in 1:length(tissueSiteDetail_)){
    suppressMessages( expTmp <- xQTLdownload_exp( genes = gene, geneType = geneType, tissueSiteDetail=tissueSiteDetail_[t], datasetId=datasetId, toSummarizedExperiment = FALSE) )
    expTmpGencodeId <- expTmp$gencodeId
    expTmp <- as.data.frame(t(expTmp))
    expTmp <- expTmp[ str_detect(rownames(expTmp), stringr::regex("^GTEX-")),, drop=FALSE]
    colnames(expTmp) <- expTmpGencodeId
    expTmp <- as.data.table(cbind(data.table(tissueSiteDetail= tissueSiteDetail_[t], sampleId = rownames(expTmp)), expTmp))
    expProfiles <- rbind(expProfiles,expTmp)
    message("== Fetching expression...",t,"/",length(tissueSiteDetail_), " - ", tissueSiteDetail_[t], " - ", nrow(expTmp)," samples.", format(Sys.time(), " | %Y-%b-%d %H:%M:%S ") )
    rm(expTmp)
  }
  message("== Done")
  expProfilesMelt <- melt( expProfiles[,-c("sampleId")], id.vars = c("tissueSiteDetail"), variable.name = "geneName", value.name = "expTPM")
  expProfilesMelt$geneName <- as.character(expProfilesMelt$geneName)
  expProfilesMelt$expTPM <- as.numeric(expProfilesMelt$expTPM)
  expProfilesMelt <- merge(expProfilesMelt, tissueSiteDetail, by="tissueSiteDetail")

  if(log10y==TRUE){
    expProfilesMelt[expTPM==0,"expTPM"]<-1
  }
  if(toTissueSite){
    tissueSiteDetail <- tissueSiteDetail[,.(tissueSite, colorHex)][,.(colorHex=colorHex[1]), by="tissueSite"]

    p1 <- ggplot(expProfilesMelt,aes(x=tissueSite, y=(expTPM), fill=tissueSite))+
      geom_violin(color="white", width=0.88, trim=FALSE, alpha=0.9, scale="width")+
      geom_boxplot(width=0.2,  alpha=0.9, outlier.size = 0.8, outlier.alpha = 0.4, outlier.shape = 21)+
      scale_fill_manual(breaks = tissueSiteDetail$tissueSite, values = tissueSiteDetail$colorHex)+
      theme_classic()+	#分组绘制
      ylab("Expression (TPM)")+
      # scale_y_log10()+
      theme(axis.text.x=element_text(size=rel(1.3), angle = 300, hjust = 0, vjust=0.5),
            axis.text.y = element_text(size=rel(1.3)),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=rel(1.3)),
            legend.position = "none"
      )
    # + geom_text(aes(x=ifelse(length(unique(genoLable$genoLabels))==3, 2, 1.5), y=max(genoLable$normExp+1.2), label=paste0("P-value: ",signif(eqtlInfo$pValue, 3)) ))
  }else{
    p1 <- ggplot(expProfilesMelt, aes(x=tissueSiteDetail, y=(expTPM), fill=tissueSiteDetail))+
      geom_violin(color="white", width=0.88, trim=FALSE, alpha=0.9, scale="width")+
      geom_boxplot(width=0.2,  alpha=0.9, outlier.size = 0.8, outlier.alpha = 0.4, outlier.shape = 21)+
      scale_fill_manual(breaks = tissueSiteDetail$tissueSiteDetail, values = tissueSiteDetail$colorHex)+
      theme_classic()+
      ylab("Expression (TPM)")+
      # scale_y_log10()+
      theme(axis.text.x=element_text(size=rel(1.3), angle = 300, hjust = 0, vjust=0.5),
            axis.text.y = element_text(size=rel(1.1)),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=rel(1.2)),
            legend.position = "none"
      )
  }
  if(log10y){
    p1 <- p1 + scale_y_log10()
  }


  print(p1)
  return(list(expProfiles=expProfiles, plot=p1))
}



#' @title Visualization of QTL specificity among multiple cells/tissues.
#'
#' @param specificityDT A data.table object from the function `xQTLanalyze_qtlSpecificity`
#' @param outPlot (character) options: "heatmap" (default) and "regression".
#' @param binNum (numeric) number of LD bins for heatmap plot. Default:4.
#' @param topTissues (numeric) number of top tissues that sorted with slope increasingly to visualize for regression plot. Default: 5
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom scales hue_pal
#' @importFrom cowplot plot_grid
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # please see function `xQTLanalyze_qtlSpecificity`
xQTLvisual_qtlSpecificity <- function(specificityDT, outPlot="heatmap", binNum=4, topTissues=5){

  # extract variable:
  snpLD <- specificityDT[['snpLD']]
  assoAllLd <- specificityDT[['assoAllLd']]
  cor_R2_logP <- specificityDT[['cor_R2_logP']]
  lm_R2_logP <- specificityDT[['lm_R2_logP']]

  data.table::setDT(snpLD)
  data.table::setDT(assoAllLd)
  data.table::setDT(cor_R2_logP)
  data.table::setDT(lm_R2_logP)


  # recut LD into bins:
  snpLD$LDbins <- as.character(cut(snpLD$R2, breaks=seq(0,1,length.out=(binNum+1)) ))
  snpLD <- snpLD[order(R2)]
  snpLD$LDorder <- 1:nrow(snpLD)
  assoAllLd <- merge(snpLD[,.(snpId=SNP_B, LDbins)], assoAllLd, by="snpId")

  # cor:
  # cor_R2_logP$corRPcut <- as.character( cut(abs(cor_R2_logP$corRP), breaks = seq(0,1,length.out=101)) )
  assoAllLd <- merge(cor_R2_logP, assoAllLd , by="tissue_label")[order(-corRP)]


  # Retain max pvalue in each bin:
  minP_f <- function(x){  data.table(tissue_label=x[1,]$tissue_label, corPR=x[1,]$corRP, LDbins= x[1,]$LDbins, logP_minMax=max(x$logP_minMax) )  }
  heatmapDT <- assoAllLd[,minP_f(.SD), by=c("LDbins", "tissue_label")]
  # fill NA bins:
  heatmapDT_allComb <- data.table::as.data.table(expand.grid(LDbins = unique(heatmapDT$LDbins), tissue_label =unique(heatmapDT$tissue_label) ))
  heatmapDT_allComb <- merge(heatmapDT_allComb, heatmapDT, by=c("LDbins", "tissue_label"), all.x = TRUE)
  heatmapDT_allComb <- merge(heatmapDT_allComb, data.table(LDbins=unique(snpLD$LDbins), ID=1:length(unique(snpLD$LDbins))), by="LDbins", all.x=TRUE)
  rm(heatmapDT)

  if(outPlot == "heatmap"){
    p1 <- ggplot(heatmapDT_allComb)+
      geom_tile(aes(x=ID, y=reorder(tissue_label, corPR), fill=logP_minMax), color="#595959")+
      scale_x_continuous(breaks = unique(heatmapDT_allComb$ID), labels = unique(heatmapDT_allComb$LDbins))+
      # geom_text(aes(x=LDorder, y=tissue_label, label = round(logP,2),  color = logP), size = 3.5)+
      # scale_fill_gradient2(low="grey", mid="orange", high="red")+
      scale_fill_gradientn(colors= colorRampPalette(c("#fcffe6", "#95de64", "#5976ba"))(length(unique(heatmapDT_allComb$logP_minMax))) )+
      theme_classic()+
      xlab("LD bins")+
      # guides(color="none")+
      theme(
        axis.text.x=element_text(size=rel(1.5), angle = 30, hjust=1, vjust=1),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(1.5)),
        axis.title.y = element_blank(),
        panel.grid.minor = element_line(colour="grey", size=0.5),
        legend.position = "top",
        legend.box="horizontal",
        legend.key.width = unit(0.06, "npc"),
        plot.title = element_text(hjust=0.5),
        plot.margin=unit(c(0.3,0,0.3,0.3),"cm")
      )+
      guides(fill = guide_colourbar(title.position="top", title.hjust = 0, title = expression(paste("Min-max normailzed  ",-log["10"],"P",sep="")) ))


     p2 <-ggplot(heatmapDT_allComb)+
      geom_tile(aes(x=1, y=reorder(tissue_label, corPR), fill=corPR),color = "black")+
      scale_x_continuous(breaks = c(1), labels = "Correlation")+
      # breaks = seq(-1,1, length.out=5), labels = seq(-1,1, length.out=5),
      scale_fill_gradientn(  colors= c("#40a9ff", "white", "#de82a7"))+
      xlab("")+
      theme_minimal()+
      theme(
        legend.position = "top",
        legend.box="horizontal",
        legend.key.width = unit(0.06, "npc"),
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=rel(1.5), angle = 30, hjust=1, vjust=1),
        panel.grid.major = element_blank(),
        plot.margin=unit(c(0.3,1,0.3,0.3),"cm")
          )+
      guides(fill = guide_colourbar(title.position="top", title.hjust = 0, title = "Correlation coefficient" ))

    p3 <- cowplot::plot_grid(p1, p2, align = "h", ncol = 2, rel_widths = c(12,1.5))
    return(p3)
  }

  # lm:
  # lm_R2_logP$colorP=colorRampPalette(c("#096dd9", "#f5f5f5", "#cf1322"))( length(unique(assoAllLd$tissue_label)) )
  lm_R2_logP$colorP <- c(  rep("#d9d9d9", nrow(lm_R2_logP)-topTissues), scales::hue_pal()(topTissues))
  regDT <- merge(lm_R2_logP, assoAllLd , by="tissue_label")[order(-slope)]
  regDT$tissue_label <- factor(regDT$tissue_label, levels = lm_R2_logP$tissue_label)
  lm_R2_logP$lineSize= 3^((seq(0.5,35,length.out=nrow(lm_R2_logP)))/10)/10

  # top tissues with max(+) and min(-) slopes for plot:
  lm_R2_logP_top <- rbind(na.omit(lm_R2_logP[slope>0][order(-slope)][1:topTissues,]) )
  lm_R2_logP_top <- cbind( lm_R2_logP_top, na.omit(rbind( na.omit(lm_R2_logP_top[slope>0][,.(x= 0.99, y=0.99*slope+intercept)]) )) )
  lm_R2_logP_top$tissue_slope <- paste0(lm_R2_logP_top$tissue_label, " (", round(lm_R2_logP_top$slope,2), ")")


  if(outPlot=="regression"){
     p <- ggplot()+
      # geom_point(aes(x=R2, y=logP_minMax, color=tissue_label))+
      geom_smooth(data=regDT, aes(x=R2, y=logP_minMax, color=tissue_label,size= tissue_label), method = "lm",se = FALSE)+
      scale_size_manual( breaks=lm_R2_logP$tissue_label, values = lm_R2_logP$lineSize )+
      scale_color_manual( breaks=lm_R2_logP$tissue_label, values =  lm_R2_logP$color )+
      scale_y_continuous( breaks=seq(0,1,0.2), labels = c("0.0",seq(0.2,0.8,0.2),"1.0") )+
      # scale_x_continuous(limits = c(0,1))+
      theme_classic()+
      ylab(expression(paste("Min-max normalized  ",-log[10],"P",sep="")))+
      xlab(expression(R^2))+
      # expand_limits(x=c(0.2, 1.8))+
      theme(
        legend.position = "none",
        axis.text = element_text(size=rel(1.4)),
        axis.title = element_text(size=rel(1.5)),
        # plot.margin = margin(0,4,0,0, "cm")
      )+
      # geom_point(data=lm_R2_logP_top, x=0.99,aes(y= intercept+slope*0.99), color="black")+
      geom_label_repel(data=lm_R2_logP_top,
                       aes(x=x, y= y, label=tissue_slope),
                       nudge_x = -0.1,
                       segment.colour="grey", segment.size = 0.5,
                       # arrow = arrow(length = unit(0.01, "npc")),
                       box.padding = 1, max.overlaps=10)

     return(p)
  }else{
    stop("\"outPlot\" can only be choosen from \"heatmap\" and \"regression\" ")
  }
}

