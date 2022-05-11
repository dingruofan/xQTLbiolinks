#' @title Plot normalized expression among genotypes for eQTL.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "snpId" or "variantId". Default: "snpId".
#' @param geneType A character string. "geneSymbol"(default) or "gencodeId".
#' @param tissueSiteDetail A character string. Can not be null.
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
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @import curl
#' @import jsonlite
#' @return A plot
#' @export
#'
#' @examples
#' \donttest{
#'  # EQTL associatons of TP53:
#'  expEqtl <- xQTLvisual_eqtlExp(variantName="rs78378222", gene ="TP53",
#'                                tissueSiteDetail="Esophagus - Mucosa")
#'  expEqtl <- xQTLvisual_eqtlExp(variantName="rs78378222", gene ="TP53",
#'                                tissueSiteDetail="Lung")
#'  expEqtl <- xQTLvisual_eqtlExp(variantName="rs3778754", gene ="IRF5",
#'                                tissueSiteDetail="Whole Blood")
#' }
xQTLvisual_eqtlExp <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8" ){
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
    message("   eQTL association was found in ", datasetId, " of gene [", gene,"] and variant [", variantName," in [", tissueSiteDetail,"].")
    message("== Done")
  }

  message("== Querying expression from API server:")
  suppressMessages(eqtlExp <- xQTLdownload_eqtlExp(variantName = eqtlInfo$snpId, gene = eqtlInfo$geneSymbol, tissueSiteDetail = tissueSiteDetail))
  if( !exists("eqtlExp") || is.null(eqtlExp) || nrow(eqtlExp)==0 ){
    stop("No expression profiles were found for gene [", gene, "] in ", tissueSiteDetail, " in ", datasetId,".")
  }else{
    message("   Normalized expression of [",nrow(eqtlExp), "] samples were obtaioned in ", tissueSiteDetail, " in ", datasetId,".")
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
    p <- ggplot(genoLable)+
      geom_boxplot( aes(x= genoLabels, y=normExp, fill=genoLabels), alpha=0.8)+
      geom_jitter(aes(x= genoLabels, y=normExp, fill=genoLabels), position=position_jitter(0.18), size=2.0, alpha=0.4, pch=21)+
      scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
      # scale_fill_manual(values=c("green", "red"))+
      # scale_fill_brewer(palette = "Dark2")+
      theme_bw()+
      labs(title = paste0(ifelse(eqtlInfo$snpId==""|| is.na(eqtlInfo$snpId), eqtlInfo$variantId, eqtlInfo$snpId), "- ", eqtlInfo$geneSymbol, " (",tissueSiteDetail,")") )+
      xlab("Genotypes")+
      ylab("Normalized expression")+
      theme(axis.text.x=element_text(size=rel(1.3)),
            axis.text.y = element_text(size=rel(1.3)),
            axis.title = element_text(size=rel(1.3)),
            legend.position = "none",
            legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid",
                                             colour ="white"),
            plot.title = element_text(hjust=0.5)
      )
    print(p)
  }

  # ggplot(genoLabelPie) +
  #   geom_bar( aes(x="", y=labelNum, fill=genoLabels), stat = "identity") + coord_polar("y", start=0)+
  #   labs(x = "", y = "", title = "") +
  #   theme_bw()+
  #   theme(
  #     axis.ticks = element_blank(),
  #     axis.text.x = element_blank(),
  #     legend.title =element_text(face="bold",size=rel(1.1)),
  #     legend.position = "right",
  #     legend.text = element_text(size=rel(1.1)),
  #     panel.border = element_blank(),
  #     panel.grid = element_blank()
  #   )+ scale_fill_discrete("geno Labels",breaks=genoLabelPie$genoLabels, label=genoLabelPie$legends)

  # gridExtra::grid.arrange()

  return(list(eqtl=eqtlInfo, exp=genoLable))
}

#' @title Plot normalized expression among genotypes for sQTL.
#'
#' @param variantName A character string. like dbsnp ID or variant id in GTEx.
#' @param gene A gene symbol or a gencode id (versioned).
#' @param variantType A character string. "snpId" or "variantId". Default: "snpId".
#' @param geneType A character string. "geneSymbol"(default) or "gencodeId".
#' @param tissueSiteDetail A character string. Can not be null.
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
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @import curl
#' @import jsonlite
#' @return Normalized expression and ggplot2 object
#' @export
#'
#' @examples
#' \donttest{
#'  # EQTL associatons of TP53:
#'  expSqtl <- xQTLvisual_sqtlExp(variantName="rs1450891501",
#'                                phenotypeId ="chr1:497299:498399:clu_54863:ENSG00000239906.1",
#'                                tissueSiteDetail="Lung")
#' }
xQTLvisual_sqtlExp <- function(variantName="", phenotypeId="", variantType="snpId", tissueSiteDetail="", datasetId="gtex_v8" ){
  genoLabels <- normExp <- labelNum <- p <- NULL

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

  # message("== Querying significant sQTL associations from API server:")
  # suppressMessages(sqtlInfo <- xQTLdownload_sqtlSig(variantName = variantName, variantType = variantType, tissueSiteDetail = tissueSiteDetail))
  # if( !exists("eqtlInfo") || is.null(eqtlInfo) || nrow(eqtlInfo)==0 ){
  #   stop("No eqtl associations were found for gene [", gene, "] and variant [", variantName,"] in ", tissueSiteDetail, " in ", datasetId,".")
  # }else{
  #   message("   eQTL association was found in ", datasetId, " of gene [", gene,"] and variant [", variantName," in [", tissueSiteDetail,"].")
  #   message("== Done")
  # }

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

  if( requireNamespace("ggplot2") ){
    p <- ggplot(genoLable)+
      geom_boxplot( aes(x= genoLabels, y=normExp, fill=genoLabels), alpha=0.8)+
      geom_jitter(aes(x= genoLabels, y=normExp, fill=genoLabels), position=position_jitter(0.18), size=2.0, alpha=0.4, pch=21)+
      scale_x_discrete( breaks=genoLableX$genoLabels, labels=genoLableX$label)+
      # scale_fill_manual(values=c("green", "red"))+
      # scale_fill_brewer(palette = "Dark2")+
      theme_bw()+
      labs(title = paste0(ifelse(varInfo$snpId==""|| is.na(varInfo$snpId), varInfo$variantId, varInfo$snpId), "- ", phenotypeId, " (",tissueSiteDetail,")") )+
      xlab("Genotypes")+
      ylab("Normalized expression")+
      theme(axis.text.x=element_text(size=rel(1.3)),
            axis.text.y = element_text(size=rel(1.3)),
            axis.title = element_text(size=rel(1.3)),
            legend.position = "none",
            legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid",
                                             colour ="white"),
            plot.title = element_text(hjust=0.5)
      )
    print(p)
  }

  return(list(varInfo=varInfo, exp=genoLable))
}



#' @title LocusZoom plot
#'
#' @param DF A data.frame or a data.table object. Four columns are required: "snps", a character list, using an rsID or chromosome coordinate (e.g. "chr7:24966446"); chromosome, chr1-chr22; Genome position; P-value.
#' @param highlightSnp Default is the snp that with lowest p-value.
#' @param population Supported population is consistent with the LDlink, which can be listed using function "LDlinkR::list_pop()"
#' @param posRange visualized genome region of interest. Default is the region that covers all snps.
#' @param token LDlink provided user token, default = NULL, register for token at https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param windowSize Window around the highlighted snp for querying linkage disequilibrium information. Default:500000
#' @param genome "grch38"(default) or "grch37".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @import utils
#' @return A data.table object and plot.
#' @export
#'
#' @examples
#' \donttest{
#'  # For GWAS:
#'  gwasFile <- tempfile(pattern = "file")
#'  gwasURL <- "https://raw.githubusercontent.com/dingruofan/exampleData/master/gwas/AD/gwasChr6Sub1.txt"
#'  utils::download.file(gwasURL, destfile=gwasFile)
#'  gwasDF <- data.table::fread(gwasFile, sep="\t", header=TRUE)
#'  gwasDF <- gwasDF[,.(rsid, chr, position,P)]
#'  xQTLvisual_locusZoom(gwasDF)
#'  xQTLvisual_locusZoom(gwasDF, posRange="chr6:3e7-7e7", population ="AFR", windowSize=200000)
#'  xQTLvisual_locusZoom(gwasDF, posRange="chr6:3e7-7e7", population ="AFR", windowSize=500000, highlightSnp="rs9271165")
#'
#'  # For eQTL:
#'  eqtlAsso <- xQTLdownload_eQTLAllAsso("RP11-385F7.1", tissueSiteDetail = "Brain - Cortex", withB37VariantId=FALSE)
#'  eqtlAsso[,c("chrom","pos")]<-rbindlist(lapply(eqtlAsso$variantId, function(x){ a=stringr::str_split(x,"_")[[1]];return(data.table(chrom=a[1], pos=a[2])) }))
#'  xQTLvisual_locusZoom( eqtlAsso[,.(snpId, chrom, pos, pValue)], population="EUR",
#'                       posRange="chr6:46488310-48376712", genomeVersion="grch38" )
#' }
xQTLvisual_locusZoom <- function( DF , highlightSnp="", population="EUR", posRange="", token="9246d2db7917", windowSize=500000, genomeVersion="grch38", snpLD=NA){
  snpId <- pos <- pValue <- logP <- pointShape<- NULL
  RS_Number <- R2 <- SNP_B <- r2Cut <- .<-NULL
  # highlightSnp=""
  # population="EUR"
  # posRange="chr6:46488310-48376712"
  # token="9246d2db7917"
  # windowSize=1e6
  # genomeVersion="grch38"

  names(DF) <- c("snpId", "chrom", "pos", "pValue")
  DF$pos <- as.integer(DF$pos)
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
    stop("Please enlarge the genomeVersion range!")
  }

  # highligth SNP:
  if(highlightSnp ==""){
    highlightSnp <- DF[which.min(pValue)]$snpId
  }

  # LD info:
  if(is.na(snpLD)){
    snpLD <- retrieveLD_ldlink(highlightSnp,population = population, windowSize = windowSize, genomeVersion = genomeVersion, token = token)
  }
  snpLD <- snpLD[,.(SNP_A=highlightSnp, SNP_B=RS_Number, R2)]

  # Set LD SNP color:
  if( nrow(snpLD)>0 ){
    # 由于多个人种存在重复LD，所以取均值：
    snpLD <- snpLD[,.(R2=mean(R2)),by=c("SNP_A","SNP_B")]
    # snpLD$colorP = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('#636363','#7fcdbb','darkgreen','#feb24c','gold'), include.lowest=TRUE))
    snpLD$r2Cut = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE))
    # snpLD$sizeP = as.character(cut(snpLD$R2,breaks=c(0,0.8, 0.9,1), labels=c(1,1.01,1.1), include.lowest=TRUE))
  }else{
    # message("No LD information of [",highlightSnp,"].")
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
                         pointColor= c("#9C8B88", "#e09351", "#df7e66", "#b75347", "#A40340"),
                         pointFill = c("#9C8B88", "#e09351", "#df7e66", "#b75347", "#096CFD"),
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
    geom_point(aes(x=pos, y=logP, fill=r2Cut, color=r2Cut, size=pointShape, shape=pointShape))+
    scale_size_manual(breaks = c('normal', "highlight"), values =  c(2,3)  )+
    scale_shape_manual(breaks = c('normal', "highlight"), values =  c(16,23) )+
    scale_color_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointColor) +
    scale_fill_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointFill) +
    # geom_text(aes(x=pos, y=logP, label=snpId ))+
    geom_label_repel(data=DF[snpId==highlightSnp,], aes(x=pos, y=logP, label=snpId) )+
    # labs(title = plotTitle )+
    xlab( xLab )+
    ylab( yLab )+
    theme_bw()+
    theme(axis.text.x=element_text(size=rel(1.3)),
          axis.title.x=element_text(size=rel(1.3)),
          axis.title.y=element_text(size=rel(1.3)),
          plot.title = element_text(hjust=0.5),
          legend.title = element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2))
    )
  if(nrow(snpLD)==0){
    p <- p+ guides( fill="none", color = "none", shape="none", size="none")
  }else{
    p <- p+ guides( shape="none", size="none", color = guide_legend(override.aes = list(size = 4)) )
  }
  print(p)
  return(p)
}

#' @title LocusCompare plot
#'
#' @param eqtlDF A data.frame or data.table with two columns: dbSNP id and p-value.
#' @param gwasDF A data.frame or data.table with two columns: dbSNP id and p-value.
#' @param highlightSnp Default is the snp that is farthest from the origin of the coordinates.
#' @param population Supported population is consistent with the LDlink, which can be listed using function LDlinkR::list_pop()
#' @param token LDlink provided user token, default = NULL, register for token at https://ldlink.nci.nih.gov/?tab=apiaccess
#' @param windowSize Window around the highlighted snp for querying linkage disequilibrium information. Default:500000
#' @param genome "grch38"(default) or "grch37".
#' @import data.table
#' @import ggplot2
#' @import stringr
#' @import ggrepel
#' @return A plot
#' @export
#'
#' @examples
#' \donttest{
#'   eqtlURL <- "https://gitee.com/stronghoney/exampleData/raw/master/eqtl/eqtlAsso.txt"
#'   gwasURL <- "https://gitee.com/stronghoney/exampleData/raw/master/gwas/AD/gwasChr6Sub1.txt"
#'   eqtlDF <- data.table::fread(rawToChar(curl::curl_fetch_memory(eqtlURL)$content), sep="\t")
#'   gwasDF <- data.table::fread(rawToChar(curl::curl_fetch_memory(gwasURL)$content), sep="\t")
#'   eqtlDF <- eqtlDF[,.(snpId, pValue)]
#'   gwasDF <- gwasDF[,.(rsid, P)]
#'   xQTLvisual_locusCompare( eqtlDF, gwasDF )
#' }
xQTLvisual_locusCompare <- function(eqtlDF, gwasDF, highlightSnp="", population="EUR",  token="9246d2db7917", windowSize=500000, genome="grch38",snpLD=NA ){
  pValue <- snpId <- distance <- logP.gwas <- logP.eqtl <- NULL
  RS_Number <- R2 <- SNP_B <- r2Cut <- pointShape<- .<-NULL
  eqtlDF <- eqtlDF[,1:2]
  gwasDF <- gwasDF[,1:2]

  setDT(eqtlDF)
  setDT(gwasDF)
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
    stop("Please enlarge the genome range!")
  }

  # log:
  DF$logP.eqtl <- (-log(DF$pValue.eqtl, 10))
  DF$logP.gwas <- (-log(DF$pValue.gwas, 10))
  DF$distance <- sqrt(DF$logP.gwas^2+DF$logP.eqtl^2)

  # highligth SNP:
  if( highlightSnp =="" ){
    highlightSnp <- DF[which.max(distance)]$snpId
  }
  # varInfo <-  xQTLquery_varId(highlightSnp)

  # get LD information:
  if(is.na(snpLD)){
    snpLD <- retrieveLD_ldlink(highlightSnp,population = population, windowSize = windowSize, genomeVersion = genome, token = token)
  }
  snpLD <- snpLD[,.(SNP_A=highlightSnp, SNP_B=RS_Number, R2)]



  # Set LD SNP color:
  if( nrow(snpLD)>0 ){
    # 由于多个人种存在重复LD，所以取均值：
    snpLD <- snpLD[,.(R2=mean(R2)),by=c("SNP_A","SNP_B")]
    # snpLD$colorP = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('#636363','#7fcdbb','darkgreen','#feb24c','gold'), include.lowest=TRUE))
    snpLD$r2Cut = as.character(cut(snpLD$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE))
    # snpLD$sizeP = as.character(cut(snpLD$R2,breaks=c(0,0.8, 0.9,1), labels=c(1,1.01,1.1), include.lowest=TRUE))
  }else{
    message("No LD information of [",highlightSnp,"].")
    snpLD <- data.table::data.table(SNP_A=character(0), SNP_B =character(0),R2=numeric(0), color=character(0),r2Cut=character(0) )
  }

  # set color:
  gwas_eqtl <- merge(DF, snpLD[,.(snpId=SNP_B, r2Cut)], by ="snpId", all.x=TRUE, sort=FALSE)
  gwas_eqtl[is.na(r2Cut),"r2Cut"]<- "(0.0-0.2]"
  gwas_eqtl[snpId==highlightSnp,"r2Cut"] <- "(0.8-1.0]"
  # set color:
  colorDT <- data.table::data.table( r2Cut = as.character(cut(c(0.2,0.4,0.6,0.8,1),breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE)),
                                     pointColor= c("#9C8B88", "#e09351", "#df7e66", "#b75347", "#A40340"),
                                     pointFill = c("#9C8B88", "#e09351", "#df7e66", "#b75347", "#096CFD"))
  colorDT <- merge(colorDT, unique(gwas_eqtl[,.(r2Cut)]), by="r2Cut",all.x=TRUE)[order(-r2Cut)]
  # gwas_eqtl <- merge( gwas_eqtl, colorDT, by="r2Cut")
  # set shape:
  gwas_eqtl$pointShape <- "normal"
  gwas_eqtl[snpId==highlightSnp,"pointShape"] <- "highlight"
  gwas_eqtl$pointShape <- as.factor(gwas_eqtl$pointShape)
  gwas_eqtl <- gwas_eqtl[order(r2Cut, logP.gwas, logP.eqtl)]
  gwas_eqtl <- rbind( gwas_eqtl[snpId!=highlightSnp ], gwas_eqtl[snpId==highlightSnp, ] )

  # title:
  plotTitle <- paste0("LocusCompare plot (",nrow(gwas_eqtl)," snps)")

  # Xlab and ylab:
  xLab <- expression(-log["10"]("p-value") (eQTL))
  yLab <- expression(-log["10"]("p-value") (GWAS))

  if( requireNamespace("ggplot2") ){
    p <- ggplot(gwas_eqtl)+
      geom_point(aes(x=logP.eqtl, y=logP.gwas, fill=r2Cut, color=r2Cut, shape=pointShape, size=pointShape))+
      scale_fill_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointFill)+
      scale_color_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointColor)+
      scale_shape_manual("Highlight",breaks = c('normal', "highlight"), values =  c(16,23) )+
      scale_size_manual("Highlight",breaks = c('normal', "highlight"), values =  c(2,3) )+
      # geom_text(aes(x=pos, y=logP, label=snpId ))+
      geom_label_repel(data=gwas_eqtl[snpId==highlightSnp,], aes(x=logP.eqtl, y=logP.gwas, label=snpId) )+
      labs(title = plotTitle )+
      xlab( xLab )+
      ylab( yLab )+
      theme_bw()+
      theme(axis.text=element_text(size=rel(1.4)),
            axis.title=element_text(size=rel(1.5)),
            plot.title = element_text(hjust=0.5),
            legend.title = element_text(size=rel(1.3)),
            legend.text = element_text(size=rel(1.2))
      )
    if( nrow(snpLD)==0){
      p <- p+ guides( fill="none", color = "none", shape="none", size="none" )
    }else{
      p <- p+ guides( shape="none", size="none", color = guide_legend(override.aes = list(size = 4)) )
    }
    print(p)
  }
  return(p)
}


#' @title Density distribution of specified genes' expression profiles in a specified tissue.
#'
#' @param genes Following gene types are supported:
#' \itemize{
#'   \item \strong{Gene symbol}.
#'
#'   A character string or a character vector (case ignored). like: "tp53","naDK","SDF4".
#'   \item \strong{Gencode/ensemble id} (versioned or unversioned).
#' }
#' @param geneType A character string. "geneSymbol"(default), "gencodeId" or "geneCategory".
#' @param tissueSiteDetail Tissue must be accessed by "tissueSiteDetailGTExv7" or "tissueSiteDetailGTExv7".
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import ggpubr
#' @importFrom SummarizedExperiment assay colData
#' @return A plot.
#' @export
#'
#' @examples
#' \donttest{
#'   genes <- c("FNDC8", "S100Z", "AQP6", "AMOT", "C3orf38", "FOXL1", "COX11", "FCN3", "DDX58", "CFI", "MS4A18", "NUDT13", "HOXA4", "VSX1")
#'   xQTLvisual_genesExp(genes, tissueSiteDetail="Lung")
#'   genes <- c("ENSG00000073598.5","ENSG00000171643.13","ENSG00000086159.12","ENSG00000126016.15","ENSG00000179021.9","ENSG00000176678.5","ENSG00000166260.10","ENSG00000142748.12","ENSG00000107201.9","ENSG00000205403.12","ENSG00000214782.7","ENSG00000166321.13","ENSG00000197576.13","ENSG00000100987.14")
#'   xQTLvisual_genesExp(genes, geneType="gencodeId", tissueSiteDetail="Liver")
#' }
xQTLvisual_genesExp <- function(genes, geneType="geneSymbol", tissueSiteDetail = "", datasetId="gtex_v8"){
  expProfiles <- xQTLdownload_exp(genes=genes, geneType = geneType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId, toSummarizedExperiment=TRUE)
  expData <- as.data.table(cbind( data.table(geneSymbol=rownames(expProfiles)), SummarizedExperiment::assay(expProfiles) ))
  expData1 <- melt(expData, id.vars="geneSymbol", variable.name = "sampleId", value.name="exp")
  ggplot( expData1, aes(x = log(exp+1,10), y = reorder(geneSymbol, -exp, median), fill = ..density..))+
    ggridges::geom_density_ridges_gradient( gradient_lwd = 1, scale = 1.4, rel_min_height = 0.05, size = 0.3) +
    # scale_fill_gradientn( colours = colorRampPalette(c("white", "blue", "red"))(27) )+
    viridis::scale_fill_viridis(name = "WCAE per SOC", option = "C")+
    xlab("Gene expression (log(TPM+1))")+
    ylab("Gene symbol")+
    theme_bw()+
    theme( legend.title = element_blank())
}


#' @title The correlation plot of two genes’ expression
#'
#' @param gene2 Gene symbol or gencode ID of two genes. Default: gene symbol.
#' @param geneType A character string. "geneSymbol"(default), "gencodeId" or "geneCategory".
#' @param groupBy Default:sex, can be choosen from pathologyNotesCategories, like: pathologyNotesCategories.mastopathy, pathologyNotesCategories.mastopathy.metaplasia.
#' @param tissueSiteDetail Tissue must be accessed by "tissueSiteDetailGTExv7" or "tissueSiteDetailGTExv7".
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#'
#' @return A plot
#' @export
#'
#' @examples
#' \donttest{
#'  gene2 = c("AMOT", "HOXA4")
#'  xQTLvisual_geneCorr( gene2, tissueSiteDetail="Liver" )
#'  xQTLvisual_geneCorr( gene2, groupBy="pathologyNotesCategories.congestion", tissueSiteDetail="Lung" )
#' }
xQTLvisual_geneCorr <- function(gene2="", geneType="geneSymbol", tissueSiteDetail = "", groupBy="sex", datasetId="gtex_v8"){
  #
  expProfiles <- xQTLdownload_exp(genes=gene2, geneType = geneType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId, toSummarizedExperiment=TRUE, pathologyNotesCategories = TRUE)
  expData <- as.data.table(cbind( data.table(geneSymbol=rownames(expProfiles)), assay(expProfiles) ))
  expData2 <- as.data.frame(t(expData[geneSymbol %in% gene2][,-c("geneSymbol")]))
  colnames(expData2) <- gene2

  sampleInfo <- as.data.table(colData(expProfiles))
  expData2 <- cbind(sampleInfo[rownames(expData2), on="sampleId"][,c("sampleId", groupBy),with=FALSE], expData2)
  expData2 <- na.omit(expData2)

  ggpubr::ggscatterhist(expData2, x=gene2[1], y=gene2[2],
                shape = 21, color = groupBy, fill=groupBy,
                margin.plot="density",
                margin.params = list(fill=groupBy, color="black", size=0.2),
                legend = c(0.9,0.15),
                ggtheme = theme_minimal())
}



#' @title eQTL significance visualization for a gene
#'
#' @param gene A gene symbol or a gencode id (versioned).
#' @param geneType A character string. "geneSymbol"(default) or "gencodeId".
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import PupillometryR
#' @return A plot
#' @export
#'
#' @examples
#' \donttest{
#'   xQTLvisual_eqtl("KIF15")
#'   xQTLvisual_eqtl("MLH1")
#' }
xQTLvisual_eqtl <- function(gene, geneType="geneSymbol", datasetId = "gtex_v8" ){
  # gene="KIF15"
  if( datasetId=="gtex_v8" ){
    gencodeVersion="v26"
  }else{
    gencodeVersion="v19"
  }

  geneInfo <- xQTLquery_gene(gene, geneType = geneType, gencodeVersion = gencodeVersion )
  geneEqtl <- xQTLdownload_eqtlSig(gene=geneInfo$geneSymbol, datasetId=datasetId)
  geneEqtlSub <- geneEqtl[,.(variantId, tissueSiteDetail, pValue)]
  geneEqtlSub$logP <- -log(geneEqtlSub$pValue, 10)
  setDF(geneEqtlSub)
  ggplot(geneEqtlSub, aes(x=reorder(tissueSiteDetail, -logP, median),y=logP))+
    PupillometryR::geom_flat_violin(data=geneEqtlSub, mapping=aes(fill=tissueSiteDetail), position=position_nudge(x=0.25), color="black", scale = "width")+
    geom_jitter(aes(color=tissueSiteDetail), width = 0.1)+
    geom_boxplot(width=0.15, position = position_nudge(x=0.25), fill="white", size=0.1)+
    coord_flip() +
    ylab(expression(-log["10"]("Pvalue")))+
    xlab("") +
    theme_bw() +
    guides(fill="none", color="none")
}


#' @title Plot distribution of gene expression among multiple tissues.
#'
#' @param genes A characer vector.
#' @param geneType "geneSymbol" or "gencodeId".
#' @param datasetId "gtex_v8" or "gtex_v7".
#' @param toTissueSite TRUE or FALSE, display all subtissues or tissue Site. Default: TURE.
#'
#' @return A data.table and a plot.
#' @export
#'
#' @examples
#' \donttest{
#'   gene = c("ENSG00000069812.11", "ENSG00000141510.16")
#'   gene="TP53"
#'   gene="HES3"
#'   geneType="gencodeId"
#'   a <- xQTLvisual_geneExpTissues("HES2",toTissueSite=TRUE)
#' }
xQTLvisual_geneExpTissues <- function(gene="", geneType="geneSymbol", datasetId="gtex_v8", toTissueSite=TRUE){
  if(datasetId == "gtex_v8"){
    tissueSiteDetail <- copy(tissueSiteDetailGTExv8)
  }else if(datasetId == "gtex_v7"){
    tissueSiteDetail <- copy(tissueSiteDetailGTExv7)
  }else{
    stop("Please choose the right datasetId!")
  }

  expProfiles <- data.table()
  tissues <- tissueSiteDetail$tissueSiteDetail
  message("== Start fetching gene expression from all tissues...")
  for(t in 1:length(tissues)){
    suppressMessages( expTmp <- xQTLdownload_exp( genes = gene, geneType = geneType, tissueSiteDetail=tissues[t], datasetId=datasetId, toSummarizedExperiment = FALSE) )
    expTmpGencodeId <- expTmp$gencodeId
    expTmp <- as.data.frame(t(expTmp))
    expTmp <- expTmp[ str_detect(rownames(expTmp), stringr::regex("^GTEX-")),, drop=FALSE]
    colnames(expTmp) <- expTmpGencodeId
    expTmp <- as.data.table(cbind(data.table(tissueSiteDetail= tissues[t], sampleId = rownames(expTmp)), expTmp))
    expProfiles <- rbind(expProfiles,expTmp)
    message("== Fetching expression ... ",t, " - ", tissues[t], " - ", nrow(expTmp)," samples." )
    rm(expTmp)
  }
  expProfilesMelt <- melt( expProfiles[,-c("sampleId")], id.vars = c("tissueSiteDetail"), variable.name = "geneName", value.name = "expTPM")
  expProfilesMelt$geneName <- as.character(expProfilesMelt$geneName)
  expProfilesMelt$expTPM <- as.numeric(expProfilesMelt$expTPM)
  expProfilesMelt <- merge(expProfilesMelt, tissueSiteDetail, by="tissueSiteDetail")

  if(toTissueSite){
    p1 <- ggplot(expProfilesMelt)+
      geom_boxplot(aes(x=reorder(tissueSite, expTPM, median), y=(expTPM), fill=tissueSite), outlier.size = 0.3)+theme_bw()+	#分组绘制
      ylab("Expression (TPM)")+
      # scale_y_log10()+
      theme(axis.text.x=element_text(size=rel(1.1), angle = 60, hjust = 1, vjust=1),
            axis.text.y = element_text(size=rel(1.1)),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=rel(1.2)),
            legend.position = "none",
            legend.background = element_rect(fill="white", size=0.5, linetype="solid",  colour ="white"),
            legend.margin = margin(0,0,0,0,unit="cm"),
            legend.title = element_blank(),
            legend.text = element_text(size=rel(1.1))
      )
  }else{
    p1 <- ggplot(expProfilesMelt)+
      geom_boxplot(aes(x=reorder(tissueSiteDetail, expTPM, median), y=(expTPM), fill=tissueSiteDetail), outlier.size = 0.3)+theme_bw()+	#分组绘制
      ylab("Expression (TPM)")+
      # scale_y_log10()+
      theme(axis.text.x=element_text(size=rel(1.1), angle = 60, hjust = 1, vjust=1),
            axis.text.y = element_text(size=rel(1.1)),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=rel(1.2)),
            legend.position = "none",
            legend.background = element_rect(fill="white", size=0.5, linetype="solid",  colour ="white"),
            legend.margin = margin(0,0,0,0,unit="cm"),
            legend.title = element_blank(),
            legend.text = element_text(size=rel(1.1))
      )
  }
  print(p1)
  return(list(expProfiles=expProfiles, plot=p1))
}
