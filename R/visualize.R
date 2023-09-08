#' @title Boxplot of normalized expression among genotypes for eQTL.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param axis_text_size (numberic) text size of the axis labels
#' @param axis_title_size (numberic) text size of the axis title
#' @param title_size (numberic) text size of the title of the plot
#' @param xlab_text Lable for x-axis
#' @param ylab_text for y-axis
#' @param ylim_v Set scale limits
#' @param title_text Title of the plot
#' @param jitter_color (A character vector) Set the point color.
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
#' \donttest{
#' expEqtl<-xQTLvisual_eqtlExp(variantName="rs3778754", gene ="IRF5",
#'                             tissueSiteDetail="Whole Blood", xlab_text="Genotypes",
#'                             ylab_text="Expression", ylim_v=c(-2,2),
#'                             axis_text_size=1.3, axis_title_size=1.3, title_size=1.4,
#'                             title_text="Genotype-expression",
#'                             jitter_color=c("#83bea5", "#e09069","#8f9dc6") )
#' }
xQTLvisual_eqtlExp <- function(variantName="", gene="", variantType="auto", geneType="auto", tissueSiteDetail="",
                               axis_text_size=1.3,axis_title_size=1.3, title_size=1.4, xlab_text="Genotypes", ylab_text="Normalized expression", ylim_v=NULL, title_text ="", jitter_color=NULL){
  genoLabels <- normExp <- labelNum <- p <- NULL
  datasetId="gtex_v8"
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
  suppressMessages(eqtlInfo <- xQTLquery_eqtl(gene = gene, variantName = variantName, geneType = geneType, variantType = variantType, tissueSiteDetail = tissueSiteDetail))
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
  # genoLabelPie <- data.table::data.table(table(genoLable$genoLabels))
  # names(genoLabelPie) <- c("genoLabels", "labelNum")
  # genoLabelPie$legends <- paste0(genoLabelPie$genoLabels, "(",genoLabelPie$labelNum,")")


  if( requireNamespace("ggplot2") ){
    if(!is.null(jitter_color)){
      jitter_color_my <- jitter_color
    }else{
      jitter_color_my <- c("#83bea5", "#e09069","#8f9dc6")
    }
    p <- ggplot( genoLable, aes(x=genoLabels, y=normExp)) +
      stat_boxplot( aes(genoLabels, normExp),
                  geom='errorbar', linetype=1, width=0.5)+  #whiskers
      geom_boxplot( aes(genoLabels, normExp),outlier.shape = NA) +
      geom_jitter(width = 0.1, aes(color=genoLabels))+ #colour="#595959",
      scale_color_manual(values = jitter_color_my[1:length(unique(genoLable$genoLabels))])+
      scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
      # labs(title = paste0(ifelse(eqtlInfo$snpId==""|| is.na(eqtlInfo$snpId), eqtlInfo$variantId, eqtlInfo$snpId), "- ", eqtlInfo$geneSymbol) )+
      labs(title = title_text)+
      xlab(xlab_text)+
      ylab(ylab_text)+
      theme_classic()+
      theme(
        # panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # axis.line = element_line(colour = "black",
        #                          linewidth = rel(1)),
        # legend.key = element_blank(),
        # strip.background = element_rect(fill = "white", colour = "black",
        #                                 linewidth = rel(2)),
        # complete = TRUE,
        axis.text.x=element_text(size=rel(axis_text_size)),
        axis.text.y = element_text(size=rel(axis_text_size)),
        axis.title = element_text(size=rel(axis_title_size)),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size = rel(title_size))
      )+
      geom_text(aes(x=ifelse(length(unique(genoLable$genoLabels))==3, 2, 1.5), y=max(genoLable$normExp+1.2), label=paste0("P-value: ",signif(eqtlInfo$pValue, 3)) ))
    if(!is.null(ylim_v)){
      p <- p+ylim(ylim_v)
    }

    # p<- ggplot( genoLable, aes(x=genoLabels, y=normExp)) +
    #   geom_violin( aes(fill=genoLabels),width=0.88, trim=FALSE, alpha=0.9, scale="width") +
    #   geom_boxplot(fill="white", width=0.2,  alpha=0.9)+
    #   scale_fill_brewer(palette="Dark2") +
    #   scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
    #   # labs(title = paste0(ifelse(eqtlInfo$snpId==""|| is.na(eqtlInfo$snpId), eqtlInfo$variantId, eqtlInfo$snpId), "- ", eqtlInfo$geneSymbol) )+
    #   xlab("Genotypes")+
    #   ylab("Normalized expression")+
    #   theme_classic()+
    #   theme(
    #     axis.text.x=element_text(size=rel(1.3)),
    #     axis.text.y = element_text(size=rel(1.3)),
    #     axis.title = element_text(size=rel(1.3)),
    #     legend.position = "none",
    #     plot.title = element_text(hjust=0.5)
    #   )+
    #   geom_text(aes(x=ifelse(length(unique(genoLable$genoLabels))==3, 2, 1.5), y=max(genoLable$normExp+1.2), label=paste0("P-value: ",signif(eqtlInfo$pValue, 3)) ))
  }
  return(list(eqtl=eqtlInfo, exp=genoLable, p=p))
}

#' @title customized Boxplot with users' own data.
#' @param genoDT (Data.framt) including two columns, "value" and "genotypes"
#' @param axis_text_size (numberic) text size of the axis labels
#' @param axis_title_size (numberic) text size of the axis title
#' @param title_size (numberic) text size of the title of the plot
#' @param xlab_text (character) Lable for x-axis
#' @param ylab_text (character) Lable for x-axis
#' @param ylim_v (numeric vector) Set scale limits
#' @param title_text (character) Title of the plot
#' @param jitter_color (A character vector) Set the point color.
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @import curl
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/eqtl/eqtlExpLabel.txt"
#' genoDT <- data.table::fread(url1)
#' box_plot <- xQTLvisual_genoBox(genoDT, title_size=1.6, title_text="Geno-Exp association" )
#' }
xQTLvisual_genoBox <- function(genoDT, axis_text_size=1.3,axis_title_size=1.3, title_size=1.4, xlab_text="Genotypes", ylab_text="Normalized expression", ylim_v=NULL, title_text ="", jitter_color=NULL){
  genoLabels <- normExp <- labelNum <- p <- NULL

  data.table::setDT(genoDT)
  names(genoDT) <- c("genotypes", "phenoValues")
  # x axis label:
  genoLableX <- data.table::as.data.table(table(genoDT$genotypes))
  names(genoLableX) <- c("genotypes", "Num")
  genoLableX$label <- paste0(genoLableX$genoLabels, "(",genoLableX$Num,")")
  genoLableX <- merge(genoDT,genoLableX, by="genotypes", sort=FALSE)

  # for Pie:
  # genoLabelPie <- data.table::data.table(table(genoLable$genoLabels))
  # names(genoLabelPie) <- c("genoLabels", "labelNum")
  # genoLabelPie$legends <- paste0(genoLabelPie$genoLabels, "(",genoLabelPie$labelNum,")")


  if( requireNamespace("ggplot2") ){
    if(!is.null(jitter_color)){
      jitter_color_my <- jitter_color
    }else{
      jitter_color_my <- c("#83bea5", "#e09069","#8f9dc6")
    }
    p <- ggplot( genoLableX, aes(x=genotypes, y=phenoValues)) +
      stat_boxplot( aes(genotypes, phenoValues),
                    geom='errorbar', linetype=1, width=0.5)+  #whiskers
      geom_boxplot( aes(genotypes, phenoValues),outlier.shape = NA) +
      geom_jitter(width = 0.1, aes(color=genotypes))+ #colour="#595959",
      scale_color_manual(values = jitter_color_my[1:length(unique(genoLableX$genotypes))])+
      scale_x_discrete( breaks=genoLableX$genotypes, labels=genoLableX$genotypes)+
      # labs(title = paste0(ifelse(eqtlInfo$snpId==""|| is.na(eqtlInfo$snpId), eqtlInfo$variantId, eqtlInfo$snpId), "- ", eqtlInfo$geneSymbol) )+
      xlab(xlab_text)+
      ylab(ylab_text)+
      labs(title = title_text)+
      theme_classic()+
      theme(
        axis.text.x=element_text(size=rel(axis_text_size)),
        axis.text.y = element_text(size=rel(axis_text_size)),
        axis.title = element_text(size=rel(axis_title_size)),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=rel(title_size))
      )
    if(!is.null(ylim_v)){
      p <- p+ylim(ylim_v)
    }
      print(p)
  }
  return(p)
}


#' @title Boxplot of normalized expression among genotypes for sQTL.
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param phenotypeId A character string. Format like: "chr1:497299:498399:clu_54863:ENSG00000239906.1"
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueSiteDetail (character) details of tissues in GTEx can be listed using `tissueSiteDetailGTExv8` or `tissueSiteDetailGTExv7`
#' @param axis_text_size (numberic) text size of the axis labels
#' @param axis_title_size (numberic) text size of the axis title
#' @param title_size (numberic) text size of the title of the plot
#' @param xlab_text (character) Lable for x-axis
#' @param ylab_text (character) Lable for x-axis
#' @param ylim_v (numeric vector) Set scale limits
#' @param title_text (character) Title of the plot
#' @param jitter_color (A character vector) Set the point color.
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
xQTLvisual_sqtlExp <- function(variantName="", phenotypeId="", variantType="auto", tissueSiteDetail="",
                               axis_text_size=1.3,axis_title_size=1.3, title_size=1.4, xlab_text="Genotypes", ylab_text="Norm.Intron-Excision Ratio", ylim_v=NULL, title_text ="", jitter_color=NULL){
  genoLabels <- normExp <-geneType<- labelNum <- p <- NULL
  datasetId="gtex_v8"
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
  varInfo <- xQTLquery_varId(variantName=variantName, variantType=variantType )


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
  # genoLabelPie <- data.table::data.table(table(genoLable$genoLabels))
  # names(genoLabelPie) <- c("genoLabels", "labelNum")
  # genoLabelPie$legends <- paste0(genoLabelPie$genoLabels, "(",genoLabelPie$labelNum,")")

  # retrieve sQTL detail:
  P_geneName <- stringr::str_split(phenotypeId, ":")[[1]][5]
  try(suppressMessages(sqtlInfo <- xQTLquery_sqtlSig(variantName = varInfo$variantId, genes = P_geneName, tissueSiteDetail=tissueSiteDetail)), silent = TRUE)
  sqtlInfo <- sqtlInfo[phenotypeId==phenotypeId]
  if(exists("sqtlInfo") & !is.null(sqtlInfo) & nrow(sqtlInfo)==1 ){
    message("Significant sQTL association found!")
    labelPvalue <- paste0("P-value: ",signif(sqtlInfo$pValue, 3))
  }else{
    labelPvalue <- ""
  }


  if( requireNamespace("ggplot2") ){

    if(!is.null(jitter_color)){
      jitter_color_my <- jitter_color
    }else{
      jitter_color_my <- c("#83bea5", "#e09069","#8f9dc6")
    }
    p <- ggplot( genoLable, aes(x=genoLabels, y=normExp)) +
      stat_boxplot( aes(genoLabels, normExp),
                    geom='errorbar', linetype=1, width=0.5)+  #whiskers
      geom_boxplot( aes(genoLabels, normExp),outlier.shape = NA) +
      geom_jitter(width = 0.1, aes(color=genoLabels))+ #colour="#595959",
      scale_color_manual(values = jitter_color_my[1:length(unique(genoLable$genoLabels))])+
      scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
      # labs(title = paste0(ifelse(eqtlInfo$snpId==""|| is.na(eqtlInfo$snpId), eqtlInfo$variantId, eqtlInfo$snpId), "- ", eqtlInfo$geneSymbol) )+
      xlab(xlab_text)+
      ylab(ylab_text)+
      labs(title = title_text)+
      theme_classic()+
      theme(
        axis.text.x=element_text(size=rel(axis_text_size)),
        axis.text.y = element_text(size=rel(axis_text_size)),
        axis.title = element_text(size=rel(axis_title_size)),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=rel(title_size))
      )+
      geom_text(aes(x=ifelse(length(unique(genoLable$genoLabels))==3, 2, 1.5), y=max(genoLable$normExp+1.2), label=paste0("P-value: ",signif(sqtlInfo$pValue, 3)) ))
    if(!is.null(ylim_v)){
      p<- p+ylim(ylim_v)
    }
    # p<- ggplot2::ggplot( genoLable, aes(x=genoLabels, y=normExp)) +
    #   geom_violin( aes(fill=genoLabels),width=0.88, trim=FALSE, alpha=0.9, scale="width") +
    #   geom_boxplot(fill="white", width=0.2,  alpha=0.9)+
    #   scale_fill_brewer(palette="Dark2") +
    #   scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
    #   # labs(title = paste0(ifelse(varInfo$snpId==""|| is.na(varInfo$snpId), varInfo$variantId, varInfo$snpId), "- ", phenotypeId, " (",tissueSiteDetail,")") )+
    #   xlab("Genotypes")+
    #   ylab("Normalized expression")+
    #   theme_classic()+
    #   theme(
    #     axis.text=element_text(size=rel(1.3)),
    #     axis.title = element_text(size=rel(1.3)),
    #     legend.position = "none",
    #     plot.title = element_text(hjust=0.5)
    #   )+
    #   geom_text(aes(x=ifelse(length(unique(genoLable$genoLabels))==3, 2, 1.5), y=max(genoLable$normExp+1.2), label= labelPvalue))
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
#' @param ylim set the minimum and maximum values for the y-axis. By default, the function will automatically determine the y-axis limits based on the data being plotted.
#' @param legend (boolean, optional) Should the legend be shown? Default: TRUE.
#' @param legend_position (character, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param point_fill (character, optional) Customized color vectors (5 kinds of colors).
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
#' \donttest{
#' library(data.table)
#' # For GWAS dataset:
#' gwasDF <- fread("https://gitee.com/stronghoney/exampleData/raw/master/gwasChr6Sub4.txt")
#' xQTLvisual_locusZoom(gwasDF)
#' # Zoom in:
#' xQTLvisual_locusZoom(gwasDF, posRange="chr6:4.7e7-4.8e7", population ="EUR")
#'
#' # For eQTL of a gene of interest (time-consuming):
#' eqtlAsso <- xQTLdownload_eqtlAllAsso("RP11-385F7.1", tissueLabel = "Brain - Cortex",
#'                                      withB37VariantId=FALSE)
#' xQTLvisual_locusZoom(eqtlAsso[,c("snpId", "chrom", "pos", "pValue")], highlightSnp="rs4711878" )
#' # Zoom in:
#' xQTLvisual_locusZoom(eqtlAsso[,c("snpId", "chrom", "pos", "pValue")], highlightSnp="rs4711878",
#'                      posRange="chr6:47.3e6-47.9e6")
#' }
xQTLvisual_locusZoom <- function( DF , highlightSnp="", population="EUR", posRange="", legend = TRUE, ylim=NULL,
                                  legend_position = c('topright','bottomright','topleft'), point_fill=NULL,
                                  snpLD=NULL){
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
    # try(snpLD <- retrieveLD_LDproxy(highlightSnp,population = population),snpLD <- snpLD[,.(SNP_A=highlightSnp, SNP_B=RS_Number, R2)])
    data.table::setDT(snpLD)

  }


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
  if(!is.null(point_fill)){
    if(length(point_fill)<5){ pointFill_my<-c(rep("blue4",5-length(point_fill)), point_fill)
    }else{ point_fill <- point_fill[1:5] }
  }else{
    pointFill_my= c('blue4','skyblue','darkgreen','orange','red')
    }
  colorDT <- data.table( r2Cut = as.character(cut(c(0.2,0.4,0.6,0.8,1),breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE)),
                         pointFill= pointFill_my,
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
  yLab <- expression(-log["10"]("p-value"))

  # xlab:
  xLab <- paste0(ifelse(stringr::str_detect(P_chrom, stringr::regex("^chr")),P_chrom, paste0("chr", P_chrom))," (",posUnit,")")


  p <- ggplot2::ggplot(DF)+
    geom_point(aes(x=pos, y=logP, fill=r2Cut,  size=pointShape, shape=pointShape), color="black")+
    scale_size_manual(breaks = c('normal', "highlight"), values =  c(3,3.5)  )+
    scale_shape_manual(breaks = c('normal', "highlight"), values =  c(21,23) )+
    # scale_y_continuous(breaks = seq(0,floor(max(DF$logP)),by=ifelse(floor(max(DF$logP))<5, 1, 5)),
    #                    labels = seq(0,floor(max(DF$logP)), by=ifelse(floor(max(DF$logP))<5, 1, 5)),limits = c(0, max(DF$logP)+0.4) )+
    # scale_color_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointColor) +
    scale_fill_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointFill) +
    # geom_text(aes(x=pos, y=logP, label=snpId ))+
    ggrepel::geom_label_repel(data=DF[snpId==highlightSnp,], aes(x=pos, y=logP, label=snpId) )+
    # labs(title = plotTitle )+
    xlab( xLab )+
    ylab( yLab )+
    theme_classic()+
    theme(axis.text=element_text(size=rel(1.3), color = "black"),
          axis.title=element_text(size=rel(1.4), color = "black"),
          plot.title = element_text(hjust=0.5),
          legend.position = "none"
    )+ guides( shape="none", color="none", size="none", fill = "none" )
  if(!is.null(ylim)){
    p <- p+ylim(ylim)
  }

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
                fill = rev(pointFill_my)) +
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
#' \donttest{
#' library(data.table)
#' # load data:
#' eqtlDF <-fread("https://gitee.com/stronghoney/exampleData/raw/master/eqtl/eqtlAsso1.txt")
#' gwasDF <-fread("https://gitee.com/stronghoney/exampleData/raw/master/gwas/AD/gwasChr6Sub3.txt")
#' # visualize:
#' xQTLvisual_locusCompare( eqtlDF, gwasDF, legend_position="topleft")
#' }
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
    # try(snpLD <- retrieveLD_LDproxy(highlightSnp,population = population), snpLD <- snpLD[,.(SNP_A=highlightSnp, SNP_B=RS_Number, R2)])
    data.table::setDT(snpLD)

  }


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
  xLab <- expression(-log["10"]("p-value") (QTL))
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
      theme(axis.text=element_text(size=rel(1.3), color = "black"),
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
#' u1 <-"http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwas/AD/gwasEqtldata.txt"
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
    # highlightSnpInfo <- xQTLquery_varId(highlightSnp)
  }else{
    message("== Highlighted SNP: [",highlightSnp,"]")
    # highlightSnpInfo <- xQTLquery_varId(highlightSnp)
  }

  # LD info:
  if( is.null(snpLD) ){
    message("== Retrieve LD information of SNP: [",highlightSnp,"]...")
    try( snpLD <- retrieveLD(DF[1,]$chrom, highlightSnp, population) )
    # try(snpLD <- retrieveLD_LDproxy(highlightSnp,population = population),snpLD <- snpLD[,.(SNP_A=highlightSnp, SNP_B=RS_Number, R2)])
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
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import viridis
#' @importFrom SummarizedExperiment assay colData
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' genes <- c("FNDC8", "S100Z", "AQP6", "AMOT", "C3orf38", "FOXL1", "COX11",
#'            "FCN3", "DDX58", "CFI", "MS4A18", "NUDT13", "HOXA4", "VSX1")
#' xQTLvisual_genesExp(genes, tissueSiteDetail="Lung")
#'
#' genes <- c("ENSG00000073598.5","ENSG00000171643.13","ENSG00000086159.12","ENSG00000126016.15",
#'            "ENSG00000179021.9","ENSG00000176678.5","ENSG00000166260.10","ENSG00000142748.12",
#'            "ENSG00000107201.9","ENSG00000205403.12","ENSG00000214782.7","ENSG00000166321.13",
#'            "ENSG00000197576.13","ENSG00000100987.14")
#' xQTLvisual_genesExp(genes, geneType="gencodeId", tissueSiteDetail="Liver")
#' }
xQTLvisual_genesExp <- function(genes, geneType="auto", tissueSiteDetail = ""){
  `..density..`<-geneSymbol <- NULL
  datasetId="gtex_v8"
  # Automatically determine the type of variable:
  if(geneType=="auto"){
    if( all(unlist(lapply(genes, function(g){ str_detect(g, "^ENSG") }))) ){
      geneType <- "gencodeId"
    }else{
      geneType <- "geneSymbol"
    }
  }

  expProfiles <- xQTLdownload_exp(genes=genes, geneType = geneType, tissueSiteDetail = tissueSiteDetail, toSummarizedExperiment=TRUE)
  expData <- as.data.table(cbind( data.table(geneSymbol=rownames(expProfiles)), SummarizedExperiment::assay(expProfiles) ))
  expData1 <- melt(expData, id.vars="geneSymbol", variable.name = "sampleId", value.name="exp")
  p <- ggplot( expData1, aes(x = log(exp+1,10), y = reorder(geneSymbol, -exp, median), fill = ..density..))+
    ggridges::geom_density_ridges_gradient( gradient_lwd = 1, scale = 1.4, rel_min_height = 0.05, size = 0.3) +
    # scale_fill_gradientn( colours = colorRampPalette(c("white", "blue", "red"))(27) )+
    viridis::scale_fill_viridis(name = "WCAE per SOC", option = "C")+
    xlab(expression("Gene expression -log"["10"]("TPM")))+
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




#' @title Box plot with jittered points for showing number and significance of eQTL associations
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' xQTLvisual_eqtl("KIF15")
#' }
xQTLvisual_eqtl <- function(gene, geneType="auto" ){
  variantId <- tissueSiteDetail <- pValue <- logP <- NULL
  . <- NULL
  datasetId = "gtex_v8"
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

  geneInfo <- xQTLquery_gene(gene, geneType = geneType )
  geneEqtl <- xQTLquery_eqtlSig(genes=geneInfo$geneSymbol)
  geneEqtlSub <- geneEqtl[,.(variantId, tissueSiteDetail, pValue)]
  geneEqtlSub$logP <- -log(geneEqtlSub$pValue, 10)
  setDF(geneEqtlSub)
  p<- ggplot(geneEqtlSub, aes(x=reorder(tissueSiteDetail, -logP, median),y=logP))+
    PupillometryR::geom_flat_violin(data=geneEqtlSub, mapping=aes(fill=tissueSiteDetail), position=position_nudge(x=0.25), color="black", scale = "width")+
    geom_jitter(aes(color=tissueSiteDetail), width = 0.1)+
    geom_boxplot(width=0.15, position = position_nudge(x=0.25), fill="white", size=0.1)+
    coord_flip() +
    ylab(expression(-log["10"]("p-value")))+
    xlab("") +
    theme_classic() +
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


#' @title Box plot with jittered points for showing number and significance of sQTL associations
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' xQTLvisual_sqtl("KIF15")
#' }
xQTLvisual_sqtl <- function(gene, geneType="auto" ){
  variantId <- tissueSiteDetail <- pValue <- logP <- NULL
  . <- NULL
  datasetId = "gtex_v8"
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

  geneInfo <- xQTLquery_gene(gene, geneType = geneType )
  geneEqtl <- xQTLquery_sqtlSig(genes=geneInfo$geneSymbol)
  geneEqtlSub <- geneEqtl[,.(variantId, tissueSiteDetail, pValue)]
  geneEqtlSub$logP <- -log(geneEqtlSub$pValue, 10)
  setDF(geneEqtlSub)
  p<- ggplot(geneEqtlSub, aes(x=reorder(tissueSiteDetail, -logP, median),y=logP))+
    PupillometryR::geom_flat_violin(data=geneEqtlSub, mapping=aes(fill=tissueSiteDetail), position=position_nudge(x=0.25), color="black", scale = "width")+
    geom_jitter(aes(color=tissueSiteDetail), width = 0.1)+
    geom_boxplot(width=0.15, position = position_nudge(x=0.25), fill="white", size=0.1)+
    coord_flip() +
    ylab(expression(-log["10"]("p-value")))+
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
#' @param log10y Display values of expression in log scale. Default: FALSE.
#' @param toTissueSite TRUE or FALSE, display all subtissues or tissue Site. Default: TURE.
#'
#' @return A list containing expression profile and a ggplot object.
#' @export
#'
#' @examples
#' \donttest{
#' # Display gene expression in specified tissues.
#' geneExpTissues <- xQTLvisual_geneExpTissues("TP53", tissues=c("Lung", "Brain","Ovary"))
#' }
xQTLvisual_geneExpTissues <- function(gene="", geneType="auto", tissues="All", log10y=FALSE, toTissueSite=FALSE){
  colorHex <- tissueSite <- expTPM <- NULL
  .<-NULL
  datasetId="gtex_v8"

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
    suppressMessages( expTmp <- xQTLdownload_exp( genes = gene, geneType = geneType, tissueSiteDetail=tissueSiteDetail_[t], toSummarizedExperiment = FALSE) )
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
            legend.position = "none",
            plot.margin=unit(c(0.3,2,0.3,0.3),"cm")
      )
  }
  if(log10y){
    p1 <- p1 + scale_y_log10()
  }


  print(p1)
  return(list(expProfiles=expProfiles, plot=p1))
}



#' @title Visualizing annotated variants
#'
#' @param snpHits A data.table object from result of xQTLanno_genomic
#' @param pValueBy Cut step of pvlaue. Defaults: 5
#' @param annoType "enrichment" or "overlapping"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwas/gwasSub.txt.gz"
#' snpInfo <- fread(url1, sep="\t")
#' snpHits <- xQTLanno_genomic(snpInfo)
#' p <- xQTLvisual_anno(snpHits)
#' }
xQTLvisual_anno <- function(snpHits, pValueBy=5, annoType="enrichment"){
  OR <- p.value <- median_OR <- logP <- cutP <- Num <- Type <- NumSum <- prop <- Freq <- snpHitsCount<- Var1 <-NULL
  .<-NULL

  if( !requireNamespace("ggforestplot") ){ stop("please install package \"ggforestplot\".") }

  typeLabel <-data.table(type=c("cpg","enhancer","promoter","exon","nonsynonymous","utr3","utr5","tfCluster", "spliceSite", "intergenic", "intron", "downstream"),
                         Type = c("CpG island", "Enhancer", "Promoter", "Exon", "Nonsynonymous", "3'UTR", "5'UTR", "TFBS", "Splice site", "Intergenic", "Intron", "Downstream"))

  if(annoType == "overlapping"){
    snpHits <- snpHits$snpHits
    snpHits$logP <- log(snpHits$pValue, base=10)*(-1)
    snpHits$logP <- ifelse(snpHits$logP==0, 1e-6, snpHits$logP)
    snpHits$cutP <- cut(snpHits$logP,breaks= ceiling(seq(0,max(snpHits$log)+pValueBy, by=pValueBy)), include.lowest=TRUE, right=FALSE )

    snpHits_count <- as.data.table(as.data.frame(table(snpHits$type)))[,.(type=Var1, Num=Freq)]
    snpHits_count <- merge(snpHits_count, typeLabel, by="type")
    snpHits_count$prop <- round(snpHits_count$Num/(sum(snpHits_count$Num))*100,2)
    snpHits_count <- snpHits_count[order(-prop)]
    snpHits_count$legendLabel=paste(snpHits_count$Type," (",snpHits_count$prop,"%",")", sep="")

    # a <- randomcoloR::distinctColorPalette(length(unique(snpHitsCount$Type)))
    # a <- randomcoloR::randomColor(length(unique(snpHitsCount$Type)))
    p1 <- ggplot(snpHits_count, aes(x = reorder(Type, Num, mean) ,y = Num,fill=Type)) +
      geom_bar(stat = "identity", position =position_dodge(0.9), width = 0.7, alpha=1)+
      # scale_y_log10( breaks=c(0.5,10^c(0:floor(log(max(snpHitsCount$Num),10)))), labels= as.character((c(0.5,10^c(0:floor(log(max(snpHitsCount$Num),10)))))) )+
      scale_fill_manual( breaks=as.character(unique(snpHits_count$Type)),
                         values=c("#D08C65","#C851DB","#DBD861","#DE80AB","#7DDDC1","#D4DDB6","#9A85D8","#8BE56B","#ABC1D8", "#e27771")[1:length(as.character(unique(snpHits_count$Type)))],
                         labels= paste0(snpHits_count$Type)
                         # labels= paste0(snpHits_count$Type, "(",snpHits_count$prop,"%)")
      )+
      theme_classic()+
      # ylim(c(0, max(snpHits_count$Num)+0.15*max(snpHits_count$Num)))+
      # geom_text(aes(label=paste(prop, "%")), position=position_dodge(width=0.9), hjust=0)+
      theme(axis.text.x = element_text(size=rel(1.3)),
            axis.text.y = element_text(size=rel(1.3)),
            axis.title = element_text(size=rel(1.3)),
            legend.title = element_text(face="bold",size=rel(1.1)),
            legend.text = element_text(size=rel(1.1)),
            legend.position = "right"
      )+
      # xlab(expression(-log["10"]("p-value")))+
      xlab("")+
      ylab("Number of variants")+coord_flip()+guides(fill="none")
    plot(p1)
    return(p1)
  }else if(annoType=="enrichment"){
    # snpHitsCount <- snpHits$snpHits[,.(num_var = nrow(.SD)),by="type"]
    # snpHitsCount$prop_var <- snpHitsCount$num_var / nrow(snpHits$gwasDF_sig)
    # snpHitsCount <- merge(snpHitsCount, typeLabel, by="type", sort=FALSE)
    snpEnrich <- snpHits$snpEnrich
    snpEnrich <- merge(snpEnrich, typeLabel, by="type", sort=FALSE)
    snpEnrich$logP <- log(snpEnrich$p.value, 10)*(-1)

    # snpEnrich_OR$Type <- factor(snpEnrich_OR[order(logP, decreasing = TRUE)]$Type )
    p1 <- ggplot(snpEnrich,aes(x=enrichment, y = reorder(Type, enrichment, median))) +
      geom_point(aes(size=logP, color=Type))+
      geom_vline(xintercept = 1, color="red", linewidth=rel(1.1), linetype="dashed", alpha=0.5)+
      ggforestplot::geom_stripes(odd = "#33333333", even = "#00000000")+
      scale_size_continuous( range = c(2,6.5))+
      theme_classic()+
      ylab("")+
      xlab("Fold Enrichment")+
      theme(axis.title = element_text(size=rel(1.2)),
            axis.text = element_text(size=rel(1.2)),
            legend.text = element_text(size=rel(1.1)))+
      guides(color="none")
    return(p1)
  }
}


#' @title Visualizing enriched variants
#'
#' @param enrichHits A data.table object from result of xQTLanno_enrich
#' @param pValueBy Cutoff of pvlaue. Defaults: 5
#' @param plotType "boxplot", or "density"
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwas/gwasSub.txt.gz"
#' snpInfo <- fread(url1, sep="\t")
#' enrichHits <- xQTLanno_enrich(snpInfo, enrichElement="TF")
#' xQTLvisual_enrich(enrichHits, plotType="density")
#' }
xQTLvisual_enrich <- function(enrichHits, pValueBy=10, plotType="boxplot"){
  cutP <- NULL
  .<-NULL

  enrichHits$logP <- log(enrichHits$pValue, base=10)*(-1)
  enrichHits$logP <- ifelse(enrichHits$logP==0, 1e-6, enrichHits$logP)
  enrichHits$cutP <- cut(enrichHits$logP,breaks= ceiling(seq(0,max(enrichHits$log)+pValueBy, by=pValueBy)), include.lowest=TRUE, right=FALSE )
  enrichHits[dist==0,"dist"]<-1

  if(plotType=="boxplot"){
    p1 <- ggplot(enrichHits)+
      geom_boxplot(aes(x=cutP,y=dist,fill=cutP))+
      scale_y_log10(breaks=10^c(0:floor(log(max(enrichHits$dist),10))), labels= as.character(as.integer(10^c(0:floor(log(max(enrichHits$dist),10))))) )+
      theme_classic()+
      xlab(expression(-log["10"]("p-value")))+
      ylab(paste("Distance",sep=""))+
      theme(axis.text.x = element_text(size = rel(1.3), angle=30, hjust=1, vjust=1),
            axis.text.y = element_text(size=rel(1.3)),
            axis.title = element_text(size=rel(1.4)),
            legend.position = "none",
            legend.title = element_blank(),
            legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid",
                                             colour ="white")
      )
    print(p1)
    return(p1)
  }else if(plotType=="density"){
    p1 <- ggplot(enrichHits, aes(x=dist))+
      geom_density(aes(color=cutP,fill=cutP),alpha=0.3)+ #添加密度图层
      # scale_fill_manual(name=expression(-log["10"]("p-value")))+
      scale_x_log10(limits = c(1,(max(enrichHits$dist)+100)),breaks=10^c(0:floor(log(max(enrichHits$dist),10))), labels= as.character(as.integer(10^c(0:floor(log(max(enrichHits$dist),10))))))+
      xlab("Distance")+
      ylab("density")+
      theme_classic()+
      theme(
        axis.text = element_text(size=rel(1.3)),
        axis.title = element_text(size=rel(1.4)),
        panel.grid = element_blank(),
        legend.position = "none"
      )
    #+ scale_x_log10()
    print(p1)
    return(p1)
  }
}

#' @title plot quantile-quantile plot with pvalue
#'
#' @param summaryDT A data.frame of one col required: pval.
#' @param legend_p TRUE or FALSE, or legend position, including: top, bottom, left and right.
#' @param binCutLogP SNPs whose logP great than this will be binned, other than not binned.
#' @param binNumber Number of bins.
#' @import ggplot2
#'
#' @export
#' @return ggplot2 object
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwas/gwasSub.txt.gz"
#' snpInfo <- data.table::fread(url1, sep="\t")
#' xQTLvisual_qqPlot(snpInfo[,.(pValue)],binCutLogP=5, binNumber=10000)
#' }
xQTLvisual_qqPlot <- function(summaryDT, legend_p=FALSE, binCutLogP=3, binNumber=1000){
  pval <- observedLogP <- expectedLogP <- NULL
  .<- NULL

  message("== Number of variants: ", nrow(summaryDT), "  ",date())
  summaryDT <- summaryDT[,1]
  names(summaryDT) <- c("pval")
  summaryDT <- na.omit(summaryDT)

  # 用原始值计算 lamdba，用bin 值作图：
  legend_p <- ifelse(legend_p, legend_p, "none")

  if(!requireNamespace("ggplot2")){stop("Require ggplot2")}

  message("calculating lamdba: ", date())
  #lamdba_p <- round(median((qnorm(summaryDT$pval/ 2))^2, na.rm =TRUE) / qchisq(0.5,1), 3)
  # http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
  lamdba_p <- median(qchisq(1-summaryDT$pval,1))/qchisq(0.5,1)

  message("== Lamdba: ", lamdba_p)

  message("== sorting... ", date())
  summaryDT <- summaryDT[order(pval)]
  message("== logging... ", date())
  summaryDT[,c("observedLogP", "expectedLogP") := .(-log10(pval), -log10(1:length(pval)/length(pval)))]

  # 要进行 bin的：
  summaryDT_cut <- summaryDT[observedLogP < binCutLogP]
  # 不进行 bin 的：
  summaryDT_noCut <- summaryDT[observedLogP >= binCutLogP]

  # 根据位置划分为 bin size 个点。
  message("== binning... ", date())
  summaryDT_cut$ID <- 1:nrow(summaryDT_cut)
  if( binNumber>=nrow(summaryDT_cut) ){ binNumber<- nrow(summaryDT_cut) }
  binSize <- round(nrow(summaryDT_cut)/binNumber,0)
  summaryDT_cut$cutF <- as.character(cut(1:nrow(summaryDT_cut), breaks=seq(0, nrow(summaryDT_cut)+binSize, binSize)))
  myFunc <- function(dt){ return(dt[c(which.max(dt$pval),which.min(dt$pval)),]) }
  myFunc2 <- function(dt){ if(nrow(dt)>5){return(dt[sample(1:nrow(dt), 5),])}else{return(dt)} }
  summaryDT_cut <-summaryDT_cut[, myFunc2(.SD), by ="cutF"]
  # summaryDT_cut <- rbind(summaryDT_cut[cutF %in% unique(summaryDT_cut$cutF)[1:10], myFunc2(.SD), by="cutF"],
  #                        summaryDT_cut[!(cutF %in% unique(summaryDT_cut$cutF)[1:10]), myFunc(.SD), by ="cutF"])

  summaryDT_all <- rbind(summaryDT_cut[,-c("ID", "cutF")], summaryDT_noCut)
  message("cutted: ",nrow(summaryDT_cut), "; non-cutted: ", nrow(summaryDT_noCut))
  rm(summaryDT_cut, summaryDT_noCut)

  message("== plotting...", date())
  P <- ggplot(summaryDT_all[,.(expectedLogP, observedLogP)])+
    geom_point(aes(x=expectedLogP, y= observedLogP ),  color="#5A90BE", size=rel(1.4), alpha=0.8)+
    # geom_line(aes(x= expectedLogP, y= observedLogP),color="#5A90BE",size=rel(1.2), alpha=0.8)+
    # scale_shape_manual(values=c(22,21,23,24))+
    geom_abline(intercept = 0, slope = 1, colour = "#f23321", size = 0.7, alpha=0.7)+
    theme_classic()+
    ylim(0, max(summaryDT_all$observedLogP)+1)+
    xlim(0, max(summaryDT_all$expectedLogP)+1)+
    theme(
      axis.text = element_text(size = rel(1.1)),
      panel.grid.major.x = element_blank(), #设置面板网格
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title=element_text(size=rel(1.3)),
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.1)),
      legend.position = legend_p,
      legend.justification = 'left',
      legend.direction='horizontal'
      # legend.margin = margin(0.1,0.1,0.1,0.1,"cm"),
      # legend.background = element_rect(fill="white", size=0.5, linetype="solid",colour ="black")
    )+xlab(expression(Expected -log["10"](p-value)))+  ylab(expression(Observed -log["10"](p-value)))+
    annotate("text", x = 0.2, y = max(summaryDT$observedLogP), label = paste0("Lamdba: ", round(lamdba_p,3)), colour = "black", size = 4.5, ,hjust=0)

  return(list(plot=P, lambda=lamdba_p))
}


#' @title Compare P-values reported in the association result file to P-values calculated from Z statistics derived from the reported beta and standard error.
#'
#' @param summaryDT A data.frame with three cols: pval,  beta, se.
#' @param binCutLogP To speed up the rendering process of the plot for tens of millions of GWAS variants, variants with a p-value below a specified threshold (binCutLogP) are randomly sampled for display.
#' @param binNumber The number of points randomly selected for plotting.
#' @param distribution_func "pnorm"(default) or "pchisq"
#' @return a list containing a data.frame of estimated pvalues and A ggplot2 object
#' @export
#'
#' @examples
#' \donttest{
#' url1 <- "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwasDFsub_MMP7.txt"
#' sumDT <- data.table::fread(url1, sep="\t")
#' xQTLvisual_PZPlot(sumDT[,.(pValue, beta, se)], distribution_func="pchisq")
#' }
xQTLvisual_PZPlot <- function(summaryDT, binCutLogP=4, binNumber=2000, distribution_func="pnorm"){
  se <- pval <- pZtest <- logPZ <- cutF <- logP <- NULL
  . <-NULL
  summaryDT <- na.omit(summaryDT)
  summaryDT <- summaryDT[,1:3]
  names(summaryDT) <- c("pval", "beta", "se")
  if(distribution_func=="pnorm"){
    summaryDT[,c("pZtest", "logP") := .(pnorm(-abs(beta/se)), log(pval, 10)*(-1))]
  }else if(distribution_func=="pchisq"){
    summaryDT[,c("pZtest", "logP") := .(pchisq((beta/se)^2, 1, lower.tail = FALSE), log(pval, 10)*(-1))]
  }

  summaryDT[,"logPZ" := log(pZtest, 10)*(-1)]
  summaryDT <- na.omit(summaryDT)

  message("== sorting... ", date())
  summaryDT <- summaryDT[order(pval)]
  message("== logging... ", date())
  # summaryDT[,c("observedLogP", "expectedLogP") := .(-log10(pval), -log10(1:length(pval)/length(pval)))]

  # 要进行 bin的：
  summaryDT_cut <- summaryDT[logP < binCutLogP]
  # 不进行 bin 的：
  summaryDT_noCut <- summaryDT[logP >= binCutLogP]

  # 根据位置划分为 bin size 个点。
  message("== binning... ", date())
  summaryDT_cut$ID <- 1:nrow(summaryDT_cut)
  if( binNumber>=nrow(summaryDT_cut) ){ binNumber<- nrow(summaryDT_cut) }
  binSize <- round(nrow(summaryDT_cut)/binNumber,0)
  summaryDT_cut$cutF <- as.character(cut(1:nrow(summaryDT_cut), breaks=seq(0, nrow(summaryDT_cut)+binSize, binSize)))
  myFunc <- function(dt){ return(dt[c(which.max(dt$pval),which.min(dt$pval)),]) }
  myFunc2 <- function(dt){ if(nrow(dt)>20){return(dt[sample(1:nrow(dt), 20),])}else{return(dt)} }
  # summaryDT_cut <-summaryDT_cut[,myFunc(.SD), by ="cutF"]
  summaryDT_cut <- rbind(summaryDT_cut[cutF %in% unique(summaryDT_cut$cutF)[1:4], myFunc2(.SD), by="cutF"],
                         summaryDT_cut[!(cutF %in% unique(summaryDT_cut$cutF)[1:4]), myFunc(.SD), by ="cutF"])

  summaryDT_all <- rbind(summaryDT_cut[,-c("ID", "cutF")], summaryDT_noCut)
  message("cutted: ",nrow(summaryDT_cut), "; non-cutted: ", nrow(summaryDT_noCut))
  rm(summaryDT_cut, summaryDT_noCut)


  p <- ggplot(summaryDT_all)+
    geom_point(aes(x=logPZ, y=logP))+
    geom_abline(intercept = 0)+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))+
    ylab(expression(Raw -log["10"](p-value)))+
    xlab(expression(Estimated -log["10"](p-value)))+
    theme_classic()+
    theme(
      axis.text = element_text(size = rel(1.3)),
      axis.title = element_text(size= rel(1.4))
    )
  print(p)
  return(list(summaryDT=summaryDT, p=p))
}


#' @title heatmap of LD-p-value of the QTL
#'
#' @param gene (character) gene symbol or gencode id (versioned or unversioned are both supported).
#' @param geneType (character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".
#' @param variantName (character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".
#' @param variantType (character) options: "auto", "snpId" or "variantId". Default: "auto".
#' @param tissueLabels (a character vector) can be listed with `ebi_study_tissues`. If is null, use all tissue / cell-types. (Default)
#' @param study (character) Studies can be listed using `ebi_study_tissues`. If is null, use all studies (Default).
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @return A list containing a data.table object and a ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' heatmapQTL <- xQTLvisual_coloc( gene="MMP7", variantName="rs11568818", study="TwinsUK")
#' }
xQTLvisual_coloc <-  function(gene="", geneType="auto", variantName="", variantType="auto", tissueLabels="", study="", population="EUR"){
  . <- NULL
  R2 <- SNP_B <- LDbins <- corRP <- ID  <- LogPpropensity <- tissue_label <- corPR <- logP_minMax <- colorRampPalette <- hue_pal <- slope <- intercept <- y <- x <- tissue_slope <- NULL
  pValue_propensity <- regDT <- lm_R2_logP_top <- label_Propensity <- pValue_eQTL <- NULL
  study_accession <- pos <- snpPanel <- variantId <- R2 <-snpId <- pValue <- tissue <- study_id <- qtl_group <- SNP_B <- logP <- NULL

  binNum <- 4
  outPlot="heatmap"

  ebi_ST <- data.table::copy(ebi_study_tissues)

  # study- tissue:
  if( all(tissueLabels!="") & study!=""){
    ebi_ST <- ebi_ST[which(tolower(tissue_label) %in% tolower(tissueLabels) & tolower(ebi_ST$study_accession) == tolower(study) ), ][,.SD[1,], by="tissue_label"]
    message(study)
  }

  # all tissue- study:
  if(study=="" & all(tissueLabels=="") ){
    ebi_ST <- ebi_ST[,.SD[1,], by="tissue_label"]
  }

  if(study!="" & all(tissueLabels=="") ){
    ebi_ST <- ebi_ST[ which(tolower(ebi_ST$study_accession) == tolower(study)),][,.SD[1,], by="tissue_label"]
  }

  if(study=="" & all(tissueLabels!="") ){
    ebi_ST <- ebi_ST[ which(tolower(tissue_label) %in% tolower(tissueLabels)),][,.SD[1,], by="tissue_label"]
  }

  if(nrow(ebi_ST)==0){
    stop("Please check study id or tissue label.")
  }

  if(gene=="" || variantName==""){
    stop("gene and variant can not be null!")
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

  # variant detail:
  variantInfo <- xQTLquery_varId(variantName, variantType = variantType)
  if( (!exists("variantInfo")) || nrow(variantInfo)==0){
    stop("Variant [",variantName,"] is not exists in GTEx, please check your input.")
  }
  if(nrow(variantInfo)>1){
    variantInfo <- variantInfo[1,]
  }

  # fetch LD:
  message("== Retrieve LD information of SNP: [",variantInfo$snpId,"]...")
  try(snpLD <- retrieveLD(variantInfo$chromosome, variantInfo$snpId, population))
  data.table::setDT(snpLD)
  snpLD <- snpLD[which(stringr::str_detect(SNP_B, stringr::regex("^rs"))), ]
  # try(snpLD <- retrieveLD_LDproxy(targetSnp= variantInfo$snpId, population = population,  windowSize = 500000,genomeVersion = "grch38", token="9246d2db7917") )
  #
  if(nrow(snpLD)<1){
    stop("No LD found for variant: [", variantInfo$snpId, "]")
  }
  message("== Number of LD-associated variants: ", nrow(snpLD))

  # 1. 获取所有组织该基因的 eQTL. 2. 为每个 LD-associated gene构建panel.
  assoAll <- data.table()
  for( i in 1:nrow(ebi_ST)){
    # 如果有多个 qtl group:
    message("")
    message("==> For tissue ",ebi_ST[i,]$tissue_label, " (",i,"/",nrow(ebi_ST),")")
    geneAsso_i <- xQTLdownload_eqtlAllAsso(gene=gene, geneType = geneType, tissueLabel = ebi_ST[i,]$tissue_label, study=ebi_ST[i,]$study_accession)
    if(!exists("geneAsso_i") ||is.null(geneAsso_i)|| nrow(geneAsso_i)==0){
      next()
    }else{
      geneAsso_i <- geneAsso_i[ qtl_group ==geneAsso_i[1,]$qtl_group,][order(pos)]
    }
    # LD-associated gene for plot:
    asso_I <- geneAsso_i[snpId %in% union(snpLD$SNP_B, snpLD$SNP_A),.(snpId, pValue, beta, tissue, tissue_label, study_id, qtl_group)]
    assoAll <- rbind(assoAll, asso_I)
    rm(asso_I)

    # all variants of gene:
    assoAllLd_i <- geneAsso_i[,.(snpId, pos, pValue)]
    assoAllLd_i <- merge(assoAllLd_i, snpLD[,.(snpId=SNP_B,R2)], by="snpId", sort=FALSE)
    assoAllLd_i$snpPanel <- paste0("s",1:nrow(assoAllLd_i))
  }

  # Retain min pvalue in each tissue_label-study, due to the duplication induced by qtl_group.
  assoAll <- assoAll[,.SD[which.min(pValue),], by=c("tissue", "tissue_label", "study_id", "qtl_group", "snpId")]
  assoAllLd <- merge(snpLD[,.(snpId=SNP_B, R2)], assoAll, by="snpId")
  assoAllLd$logP <- ifelse(assoAllLd$pValue==0, 0, (-log(assoAllLd$pValue,10)))
  # scale logP by tissue group (min-max method):
  assoAllLd <- assoAllLd[,.(snpId, R2, beta, logP, logP_minMax=(logP-min(logP))/(max(logP)-min(logP))*2-1 ),by="tissue_label"]

  # cor:
  cor_R2_logP <- assoAllLd[,.(corRP=cor(R2, logP_minMax), corPvalue=cor.test(R2, logP_minMax)$p.value),by="tissue_label"][order(corRP)]
  cor_R2_logP$logCorP <- log10( cor_R2_logP$corPvalue)*(-1)

  # extract variable:
  # snpLD <- propensityRes[['snpLD']]
  # assoAllLd <- propensityRes[['assoAllLd']]
  # cor_R2_logP <- propensityRes[['cor_R2_logP']]


  data.table::setDT(snpLD)
  data.table::setDT(assoAllLd)
  data.table::setDT(cor_R2_logP)

  # topTissues <- nrow(cor_R2_logP[corRP <= cor_cutofff])

  # recut LD into bins:
  snpLD$LDbins <- as.character(cut(snpLD$R2, breaks=seq(0,1,length.out=(binNum+1)) ))
  snpLD <- snpLD[order(R2)]
  snpLD$LDorder <- 1:nrow(snpLD)
  assoAllLd <- merge(snpLD[,.(snpId=SNP_B, LDbins)], assoAllLd, by="snpId")

  # cor:
  # cor_R2_logP$corRPcut <- as.character( cut(abs(cor_R2_logP$corRP), breaks = seq(0,1,length.out=101)) )
  assoAllLd <- merge(cor_R2_logP, assoAllLd , by="tissue_label")[order(-corRP)]


  # Retain max pvalue in each bin:
  minP_f <- function(x){  data.table::data.table(tissue_label=x[1,]$tissue_label, corPR=x[1,]$corRP, logCorP=x[1,]$logCorP, LDbins= x[1,]$LDbins, logP_minMax=max(x$logP_minMax) )  }
  heatmapDT <- assoAllLd[,minP_f(.SD), by=c("LDbins", "tissue_label")]
  # fill NA bins:
  heatmapDT_allComb <- data.table::as.data.table(expand.grid(LDbins = unique(heatmapDT$LDbins), tissue_label =unique(heatmapDT$tissue_label) ))
  heatmapDT_allComb <- merge(heatmapDT_allComb, heatmapDT, by=c("LDbins", "tissue_label"), all.x = TRUE)
  heatmapDT_allComb <- merge(heatmapDT_allComb, data.table(LDbins=unique(snpLD$LDbins), ID=1:length(unique(snpLD$LDbins))), by="LDbins", all.x=TRUE)
  if(nrow(heatmapDT_allComb)==0){
    stop("Not enough data for plot!")
  }
  rm(heatmapDT)
  heatmapDT_allComb$tissue_label <- factor(heatmapDT_allComb$tissue_label, levels = unique(heatmapDT_allComb[order(corPR)]$tissue_label))

  if(outPlot == "heatmap"){
    p1 <- ggplot(heatmapDT_allComb)+
      geom_tile(aes(x=ID, y=tissue_label, fill=logP_minMax), color="#595959")+
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
        # panel.spacing.x = unit(15, 'mm'),
        legend.direction = "horizontal",
        legend.key.width = unit(0.065, "npc"),
        plot.title = element_text(hjust=0.5),
        plot.margin=unit(c(2.5,0,0.3,0.3),"cm")
      )+
      guides(fill = guide_colourbar(title.position="top", title.hjust = 0, title = expression(paste("Normailzed  ",-log["10"],p-value,sep="")) ))


    p2 <- ggplot(unique(heatmapDT_allComb[,.(tissue_label, corPR)]))+
      geom_bar(aes(x=corPR, y=tissue_label, fill=corPR), color="black", stat="identity")+
      # geom_tile(aes(x=1, y=reorder(tissue_label, corPR), fill=corPR),color = "black")+
      # scale_x_continuous(breaks = c(0.5), labels = "Correlation")+
      # geom_text(aes(x=1, y=reorder(tissue_label, LogPpropensity),label = round(heatmapDT_allComb$LogPpropensity,2)), color="#595959")+
      # breaks = seq(-1,1, length.out=5), labels = seq(-1,1, length.out=5),
      # scale_fill_gradientn(colors= c("#66bfdf", "white", "#ea5a5a") )+
      scale_fill_gradient2(low= "#03b1f0", mid="#f6ffed", high="#ea5a5a", midpoint = 0)+
      xlab("")+
      theme_classic()+
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(0.065, "npc"),
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=rel(1.5), angle =90, hjust=1, vjust=1),
        panel.grid.major = element_blank(),
        plot.margin=unit(c(2.5,1,0.3,0.3),"cm")
      )+
      guides(fill = guide_colourbar(title.position="top", title.hjust = 0, title = "Pearson Correlation" ))

    p3 <- cowplot::plot_grid(p1, p2, align = "h", ncol = 3, rel_widths = c(12,3,1))
    print(p3)
    return(list(assoAllLd=assoAllLd,p=p3))
  }
}




#' @title Advanced CM plot of 11
#' @param gwasDF A data.frame of summary statistics data, including four cols arranged in the following order: SNP name, chomosome, position, p-value.
#' @param pvalue_cutoff Default: 1e-4. The Manhattan plot is a helpful tool for visualizing genome-wide association study results. However, when there are a large number of SNPs, the plot can become difficult to render and generate a large file size. This is due to the stacking of non-significant SNPs at the bottom of the plot. To address this issue, we can choose to filter out some of the non-significant SNPs or randomly select a subset of them to plot. This will improve the readability of the plot and reduce the file size.
#' @param num_snp_selected Default: 2000. Number of SNPs randomly selected for each chromosome.
#' @return A pdf format figure.
#' @export
#'
#' @examples
#' \donttest{
#' gwasDF <- data.table::fread(
#'       "http://bioinfo.szbl.ac.cn/xQTL_biolinks/xqtl_data/gwas/AD/gwasChr6Sub.txt")
#' xQTLvisual_manhattan(gwasDF[,.(rsid, chr, position,P)])
#' }
xQTLvisual_manhattan <- function(gwasDF, pvalue_cutoff=1e-4, num_snp_selected=2000){
  pValue <- Chromosome <-
  if( !requireNamespace("CMplot") ){ stop("please install package \"CMplot\".") }
  data.table::setDT(gwasDF)
  gwasDF <- gwasDF[,1:4]
  names(gwasDF) <- c("SNP", "Chromosome", "pos", "pValue")
  chroms<- unique(gwasDF$Chromosome)
  gwasDF_CMp_1 <- gwasDF[ pValue < pvalue_cutoff]
  gwasDF_CMp_2 <- data.table()
  for(i in 1:length(chroms)){
    message(i)
    a <- gwasDF[ Chromosome == chroms[i] & pValue>=pvalue_cutoff]
    if(nrow(a)>num_snp_selected){
      gwasDF_CMp_2 <- rbind(gwasDF_CMp_2, a[sample(1:nrow(a), num_snp_selected)])
    }else{
      gwasDF_CMp_2 <- rbind(gwasDF_CMp_2, a)
    }
  }
  gwasDF_CMp<- rbind(gwasDF_CMp_1, gwasDF_CMp_2)

  minP_eachChr <- gwasDF_CMp[,.SD[which.min(pValue)],by="Chromosome"]
  minP_eachChr$snp_color <- "red"

  CMplot::CMplot(gwasDF_CMp, plot.type="m",LOG10=TRUE,
         col=c("grey30","grey60"),
         # col=c("#4197d8", "#f8c120", "#413496", "#495226","#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"),
         highlight=minP_eachChr$SNP,
         highlight.col=minP_eachChr$snp_color,
         highlight.text=minP_eachChr$SNP,
         highlight.cex=1.7,highlight.pch=c(17),
         highlight.text.col="black", threshold=5e-08, threshold.lty=2,  signal.cex=0.7, #file.name="GWAS_breast_cancer/CG0050/Manhattan.pValue.pdf",
         amplify=FALSE,file="pdf",dpi=600,file.output=TRUE,verbose=TRUE,width=14,height=6)

}


