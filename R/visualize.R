#' @title Plot normalize expression among genotypes.
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
#'
#'
#'   dot num+++
#'  # EQTL associatons of TP53:
#'  GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Esophagus - Mucosa")
#'  GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Lung")
#'  GTExvisual_eqtlExp(variantName="rs3778754", gene ="IRF5", tissueSiteDetail="Whole Blood")
#' }
#'
#'
GTExvisual_eqtlExp <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8" ){
  genoLabels <- normExp <- labelNum <- p <- NULL
  # variantName="rs78378222"
  # gene="TP53"
  # tissueSiteDetail="Esophagus - Mucosa"
  # datasetId="gtex_v8"
  # variantType="snpId"
  # geneType="geneSymbol"

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
  suppressMessages(eqtlInfo <- GTExdownload_eqtlAll(gene = gene, variantName = variantName, geneType = geneType, variantType = variantType, tissueSiteDetail = tissueSiteDetail))
  if( !exists("eqtlInfo") || is.null(eqtlInfo) || nrow(eqtlInfo)==0 ){
    stop("No eqtl associations were found for gene [", gene, "] and variant [", variantName,"] in ", tissueSiteDetail, " in ", datasetId,".")
  }else{
    message("   eQTL association was found in ", datasetId, " of gene [", gene,"] and variant [", variantName," in [", tissueSiteDetail,"].")
    message("== Done")
  }

  message("== Querying expression from API server:")
  suppressMessages(eqtlExp <- GTExdownload_eqtlExp(variantName = eqtlInfo$snpId, gene = eqtlInfo$geneSymbol, tissueSiteDetail = tissueSiteDetail))
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

  stopifnot(require(ggplot2))
  stopifnot(require(ggrepel))
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

  print(p)
  return(list(eqtl=eqtlInfo, exp=genoLable))
}


#  注意，这里使用了 dbsnp api 来获取一定基因组范围内的 snp，GTEx v8 用的是 dbsnp 151， GTEx v7 用的是 dbsnp 147. 详见：https://www.gtexportal.org/home/datasets; https://zenodo.org/record/3572799/files/GTEx-V7_HapMap-2017-11-29_README.txt
# 用户设定范围：1. SNP 上下游，2. gene 上下游 3. 自定义范围。
# 根据范围获得GWAS的SNP。
# 根据GWAS SNP获得 eqtl 信息。
#' @title significance analysis of GTEX and eQTL  NO.
#'
#' @param gwasDF A data.frame containing GWAS summary info. At least four columns are expected, including "snpId", "chrom", "pos" and "pValue".
#' @param queryTerm Term of interest, which can be a gene (like "ABCB9", "ENSG00000141510.16"), a variant (like "rs7953894", "chr12_122920419_C_A_b38") or a genome coordinate (default: "chr1:1-300000").
#' @param queryType A character string that stands for the type of \"queryTerm\". For gene, "geneSymbol" or "gencodeId", for genome coordinate, "coordinate".
#' @param eqtlTrait eqtl trait. If this parameter is not specified, it is the same as the "queryTerm".
#' @param eqtlType eqtl type. If this parameter is not specified, it is the same as the "queryType".
#' @param tissueSiteDetail tissue
#' @param flankUp A whole integer. Upstream flanking distance of queried term. Uint: KB. Default: 100. Note: when "queryType" is "coordinate", this parameter is ignored.
#' @param flankDown A whole integer. Downstream flanking distance of queried term. Uint: KB. Default: 100. Note: when "queryType" is "coordinate", this parameter is ignored.
#' @param datasetId A character string. "gtex_v8" or "gtex_v7". Default: "gtex_v8".
#' @import data.table
#' @import curl
#' @import stringr
#' @import jsonlite
#' @import ggplot2
#' @import ggrepel
#' @return A plot.
#' @export
#'
#' @examples
#' \donttest{
#'  gwasDF_raw <- data.table::fread("D:\\R_project\\GWAS_Type-2-Diabetes_Wood_2016.txt.gz", sep="\t", header=TRUE)
#'  gwasDF <- gwasDF_raw[,.(rsid, pvalue)]
#'
#'
#'  GTExanalyze_eqtlGWAS(gwasDF, traitGene="SFTPD", traitGeneType="geneSymbol",
#'                       tissueSiteDetail="Artery - Tibial", study_id="gtex_v8")
#'  GTExanalyze_eqtlGWAS(gwasDF, traitGene="rs12778583", traitGeneType="snpId",
#'                       tissueSiteDetail="Thyroid", study_id="gtex_v7")
#'  GTExanalyze_eqtlGWAS(gwasDF, traitGene="12:122404966-124404966", traitGeneType="coordinate",
#'                       tissueSiteDetail="Whole Blood", study_id="gtex_v7")
#'
#'  traitGene="SFTPD"
#'  tissueSiteDetail="Artery - Tibial"
#'  geneAsso<- GTExdownload_assoAll(traitGene, tissueSiteDetail=tissueSiteDetail)
#'
#'  gwasDF <- data.table::fread("D:\\R_project\\GWAS_Type-2-Diabetes_Wood_2016.txt.gz", sep="\t", header=TRUE)
#'  gwasDF <- gwasDF[,.(rsid, pvalue, Freq)]
#'  colocResult <- coloc.abf(dataset1=list(pvalues=gwas_eqtl$pValue.gwas, type="cc", S=0.33, N=nrow(gwasDF)),
#'  dataset2=list(pvalues=gwas_eqtl$pValue.etql, type="quant", N=nrow(eqtlAsso)), MAF=gwas_eqtl$maf)
#'
#' }
GTExvisual_eqtlGWAS <- function(gwasDF, traitGene="", traitGeneType="geneSymbol", tissueSiteDetail="", study_id="gtex_v8"){
  pos <- pValue.gwas <- pValue.etql <-NULL
  # gwasDF <- data.table::fread("../GWAS_White-Blood-Cell-Traits_Tajuddin_2016.txt", sep="\t", header=TRUE)
  # gwasDF <- data.table::fread("../GWAS_Type-2-Diabetes_Wood_2016.txt.gz", sep="\t", header=TRUE)
  # gwasDF <- gwasDF[trait=="White-Blood-Cell-Counts",][,.(rsid, chr, POS)]
  # gwasDF <- gwasDF[,.(rsid, chr, snp_pos, pvalue, effect_allele, non_effect_allele)]

  # traitGene="ABCB9"
  # traitGeneType="geneSymbol"
  # tissueSiteDetail="Whole Blood"
  # study_id="gtex_v8"

  # queryTerm="rs7953894"
  # queryType ="snpId"

  cutNum <- 100
  gencodeVersion <- "v26"
  # rename:
  names(gwasDF) <- c("snpId", "pValue")

  # parameter check: snpId
  message("==Checking parameter!")
  # random check:
  if( nrow(gwasDF)>10000){
    set.seed(1)
    randomId <- sample(nrow(gwasDF),10000,replace = FALSE)
    if( !all(unlist(lapply(gwasDF$snpId[randomId], function(x){ stringr::str_detect(x,stringr::regex("^rs[0-9]{1,30}[0-9]$")) }))) ){
      message("The first column of \"gwasDF\" must be snp ID, which start with \"rs\"! ")
      return(data.table::data.table())
    }
    rm(randomId)
  }
  # parameter check: gwasDF : chromesome
  # if( study_id=="gtex_v8"){
  #   gencodeVersion <- "v26"
  #   genomeBuild="GRCh38/hg38"
  #   if( !all( unlist( lapply( unique(gwasDF$chrom), function(x){ stringr::str_detect(x, stringr::regex("^chr1[0-9]$|^chr2[0-3]$|^chr[1-9]$|^chr[xXyY]$")) })) ) ){
  #     gwasDF$chrom <- paste0("chr",gwasDF$chrom)
  #     message("For gtex_v8, the second column of \"gwasDF\" is chromosome, which can chosen from the \"chr1-chr23, chrX, chrY\"."," \"chr\" is added as prefix. ")
  #   }
  # }else{
  #   message("Parameter \"study_id\" must be \"gtex_v8\"  ")
  #   return(data.table::data.table())
  # }
  # parameter check: pos
  # if( any(!is.wholenumber(gwasDF$pos)) || any(gwasDF$pos<=0) ){
  #   stop("The third column of \"gwasDF\" must be an positive integer!")
  # }
  # parameter check: traitGeneType
  if( is.null(traitGeneType) ||  any(is.na(traitGeneType)) || any(traitGeneType=="") || length(traitGeneType)!=1 ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }else if( !(traitGeneType %in% c("geneSymbol", "gencodeId")) ){
    stop("Parameter \"geneType\" should be choosen from \"geneSymbol\", \"gencodeId\".")
  }
  message("== Done")

  #
  message("== Querying gene info from API server:")
  # check network:
  bestFetchMethod <- apiAdmin_ping()
  if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
    # message("Note: API server is busy or your network has latency, please try again later.")
    return(NULL)
  }
  suppressMessages( geneInfo <- GTExquery_gene(traitGene, geneType = traitGeneType, gencodeVersion = gencodeVersion) )
  if( exists("geneInfo") && nrow(geneInfo)>0 ){
    message("== Done")
  }else if( nrow(geneInfo)>1 ){
    message("Overall ",nrow(geneInfo)," genes were detected: [",paste0(geneInfo$gencodeId, collapse = ","),"]. ", "Please select one of the genes and rerun the function.")
  }else{
    message("Can not get gene info for [", traitGene,"] in ", study_id,".")
  }

  ######## Fetch all associations:
  eqtlAsso <- GTExdownload_assoAll(gene=geneInfo$gencodeId, geneType = "gencodeId", tissueSiteDetail = tissueSiteDetail )
  # if no eqtl got:
  if( nrow(eqtlAsso)==0 || !exists("eqtlAsso")){
    stop("No eqtl associations were found for ",traitGeneType,": [",traitGene,"].")
  }

  # merge and obtain info：
  gwas_eqtl <- merge(gwasDF, eqtlAsso, by="snpId", sort=FALSE, suffixes = c(".gwas",".etql"))
  if( nrow(gwas_eqtl)==0 || !exists("gwas_eqtl") ){
    message("No intersection of GWAS and eQTL dataset!")
    return(data.table::data.table())
  }

  result <- coloc.abf(dataset1=list(pvalues=gwas_eqtl$pValue.gwas, type="cc", s=0.33, N=50000),
                      dataset2=list(pvalues=gwas_eqtl$pValue.etql, type="quant", N=10000), MAF=gwas_eqtl$maf)


  ######## Fetch all LD info:
  geneLD <- GTExdownload_ld(gene=geneInfo$gencodeId, geneType = "gencodeId", datasetId = study_id, recordPerChunk = recordPerChunk)


  gwas_eqtlLD <- gwas_eqtl[,.(snpId, pValue.gwas, pValue.etql,variant)]

  if( requireNamespace("ggplot2") ){
    p<- ggplot(gwas_eqtl[,.(pValue.gwas, pValue.etql)])+
      geom_point( aes(x= (-1) * log(pValue.gwas,10), y = (-1)* log(pValue.etql,10)) )+
      theme_bw()+
      theme(axis.title.x=element_text(size=rel(1.3)),
            axis.title.y=element_text(size=rel(1.3)),
            legend.title = element_text(face="bold",size=rel(1.1)),
            legend.text = element_text(size=rel(1.1)),
            legend.position = "none",
      )
    print(p)
  }

  return(p)

}


#' @title eql trait plot
#'
#' @param gene gene name/gencode ID
#' @param geneType "geneSymbol" or "gencodeId"
#' @param tissueSiteDetail gene name
#' @param datasetId "gtex_v7" or "gtex_v8"
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @return A data.table
#' @export
#'
#' @examples
#' \donttest{
#'  GTExvisual_eqtlTrait("ENSG00000112137.12", "gencodeId",tissueSiteDetail = "Adipose - Subcutaneous", datasetId="gtex_v7")
#'  GTExvisual_eqtlTrait("tp53",tissueSiteDetail = "Adipose - Subcutaneous", datasetId="gtex_v7")
#' }
GTExvisual_eqtlTrait <- function(gene="", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8"){

  # gene = "ENSG00000112137.12"
  # geneType = "gencodeId"
  # tissueSiteDetail = "Adipose - Subcutaneous"
  # datasetId = "gtex_v7"

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
  suppressMessages( eqtlOfgene <- GTExdownload_eqtlAll( gene=gene, geneType = geneType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId))
  if( !exists("eqtlOfgene") || is.null(eqtlOfgene) || nrow(eqtlOfgene)==0 ){
    stop("No eqtl associations were found for gene [", gene, "] in ", tissueSiteDetail, " in ", datasetId,".")
  }else{
    message("   Totally, ", nrow(eqtlOfgene)," eQTL association was found in ", datasetId, " of gene [", gene,"] in [", tissueSiteDetail,"].")
    message("== Done")
  }

  if( nrow(eqtlOfgene)<=2 ){
    warning("Only ",nrow(eqtlOfgene)," eqtl ", ifelse(nrow(eqtlOfgene)==1,"association was","associations were")," detected, please extend the genome range using parameter \"coordinate\". " )
    return(NULL)
  }
  # split coordinate from variantId:
  eqtlOfgene <- cbind(eqtlOfgene, data.table::rbindlist( lapply(eqtlOfgene$variantId, function(x){ splitOut =stringr::str_split(x, stringr::fixed("_"))[[1]]; data.table(chrom=splitOut[1], pos=splitOut[2]) } )  ))
  eqtlOfgene$logP <- (-1*log(eqtlOfgene$pValue, 10))
  eqtlOfgene$pos <- as.integer(eqtlOfgene$pos)
  # order by logP desc:
  eqtlOfgene <- eqtlOfgene[order(-logP)]

  # title:
  plotTitle <- paste0(gene, " (",ifelse(datasetId=="gtex_v8",unique(eqtlOfgene$chrom), paste0("chr",unique(eqtlOfgene$chrom))),
                      ":",
                      paste0(range(eqtlOfgene$pos), collapse = "-")
                      ,")")

  # ylab and unit:
  posUnit <- "Bb"
  if( any(range(eqtlOfgene$pos)>10^6)){
    eqtlOfgene$pos <- eqtlOfgene$pos/10^6
    posUnit <- "Mb"
  }else if( all(range(eqtlOfgene$pos)<10^6) && all(range(eqtlOfgene$pos)>10^3) ){
    eqtlOfgene$pos <- eqtlOfgene$pos/10^3
    posUnit <- "Kb"
  }else{
    posUnit <- "Bb"
  }
  yLab <-expression(-log["10"]("Pvalue"))

  # xlab:
  xLab <- paste0(ifelse(datasetId=="gtex_v8",unique(eqtlOfgene$chrom),paste0("chr",unique(eqtlOfgene$chrom)))," (",posUnit,")")

  # color: 高亮最显著的那一个，并标注。
  eqtlOfgene$colorP <- "grey"
  eqtlOfgene$sizeP <- 2
  eqtlOfgene[which(eqtlOfgene$logP>=quantile(eqtlOfgene$logP, seq(0,1,0.1))[10] & eqtlOfgene$logP<quantile(eqtlOfgene$logP, seq(0,1,0.1))[11] ),]$colorP <- "green"
  eqtlOfgene[which(eqtlOfgene$logP>=quantile(eqtlOfgene$logP, seq(0,1,0.1))[10] & eqtlOfgene$logP<quantile(eqtlOfgene$logP, seq(0,1,0.1))[11] ),]$sizeP <- 2
  eqtlOfgene[1,]$colorP <- "red"
  eqtlOfgene[1,]$sizeP <- 2

  stopifnot(require(ggplot2))
  stopifnot(require(ggrepel))
  p <- ggplot(eqtlOfgene)+
    geom_point(aes(x=pos, y=logP, color=colorP, size=sizeP))+
    scale_color_manual(breaks = eqtlOfgene$colorP, values = eqtlOfgene$colorP)+
    # scale_size_manual(breaks = as.character(eqtlOfgene$sizeP), values =  as.numeric(eqtlOfgene$sizeP) )+
    # geom_text(aes(x=pos, y=logP, label=snpId ))+
    geom_label_repel(data=eqtlOfgene[1:1,], aes(x=pos, y=logP, label=snpId) )+
    labs(title = plotTitle )+
    xlab( xLab )+
    ylab( yLab )+
    theme_bw()+
    theme(axis.text.x=element_text(size=rel(1.3)),
          axis.title.x=element_text(size=rel(1.3)),
          axis.title.y=element_text(size=rel(1.3)),
          plot.title = element_text(hjust=0.5)
    )+guides(color="none", size="none")
  print(p)
  p


}






