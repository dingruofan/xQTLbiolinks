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
#'
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
#' @param traitGene eqtl trait. If this parameter is not specified, it is the same as the "queryTerm".
#' @param traitGeneType eqtl type. If this parameter is not specified, it is the same as the "queryType".
#' @param tissueSiteDetail tissue
#' @param study_id A character string. Default: "gtex_v8".
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
#'
#'
#' }
GTExvisual_eqtlGWAS <- function(gwasDF, traitGene="", traitGeneType="geneSymbol", tissueSiteDetail="", study_id="gtex_v8"){
  .<-NULL
  snpId <- variant <- logP <- pos <- colorP <- sizeP <- snpId <- NULL
  pos <- pValue.gwas <- pValue.etql <-NULL
  # gwasDF <- data.table::fread("../GWAS_White-Blood-Cell-Traits_Tajuddin_2016.txt", sep="\t", header=TRUE)
  # gwasDF <- data.table::fread("../GWAS_Type-2-Diabetes_Wood_2016.txt.gz", sep="\t", header=TRUE)
  # gwasDF <- gwasDF[trait=="White-Blood-Cell-Counts",][,.(rsid, chr, POS)]
  # gwasDF <- gwasDF[,.(rsid, chr, snp_pos, pvalue, effect_allele, non_effect_allele)]

  # traitGene="LPAR2"
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
  ######## Fetch all LD info:
  snpSig <- "rs11671619"
  snpEqtl <- GTExdownload_eqtlSig(snpSig)
  geneLD <- data.table::data.table()
  for(i in unique(snpEqtl$gencodeId)){
    message(i)
    suppressMessages( geneLDTmp <- GTExdownload_ld(gene=i, geneType = "gencodeId", datasetId = study_id, recordPerChunk = recordPerChunk) )
    geneLD <- rbind(geneLD, geneLDTmp)
  }
  snpSigLD <- rbind(geneLD[snpId_1==snpSig,], geneLD[snpId_2==snpSig,])

  geneLD <- GTExdownload_ld(gene=geneInfo$gencodeId, geneType = "gencodeId", datasetId = study_id, recordPerChunk = recordPerChunk)
  snpSig <- "rs11671619"
  gwas_eqtl <- merge(gwas_eqtl, )

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
#'  GTExvisual_eqtlTrait("ENSG00000112137.12", "gencodeId",
#'    tissueSiteDetail = "Adipose - Subcutaneous", datasetId="gtex_v7")
#'  GTExvisual_eqtlTrait("tp53",
#'    tissueSiteDetail = "Adipose - Subcutaneous", datasetId="gtex_v7")
#' }
GTExvisual_eqtlTrait <- function(gene="", geneType="geneSymbol", highlightSnp="", tissueSiteDetail="", study ="gtex_v8"){
  .<-NULL
  logP <- pos <- colorP <- sizeP <- snpId<-NULL
  # gene = "ENSG00000112137"
  # geneType = "gencodeId"
  # tissueSiteDetail = "Adipose - Subcutaneous"
  # study = "gtex_v8"

  gencodeVersion <- "v26"
  if( study == "gtex_v8" ){
    gencodeVersion <- "v26"
    tissueSiteDetailGTEx <- data.table::copy(tissueSiteDetailGTExv8)
  }

  # check gene:
  if( length(gene) >1 || !is.character(gene) ){
    stop("Parameter \"gene\" must be a character string. ")
  }else if(gene ==""){
    stop("Parameter \"gene\" can not be null. ")
  }

  # check tissueSiteDetail:
  # import tissueSiteDetailGTEx according to datasetId
  qtl_groups <- EBIquery_allTerm("qtl_groups")
  qtl_tissue <- merge( tissueSiteDetailGTEx,qtl_groups, by.x="tissueSiteDetailId", by.y="qtl_group")
  setDT(qtl_tissue)
  # check tissueSiteDetail:
  if( is.null(tissueSiteDetail) ||  any(is.na(tissueSiteDetail)) || tissueSiteDetail==""   ){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if(length(tissueSiteDetail)!=1){
    stop("Parameter \"tissueSiteDetail\" should be chosen from following tissue names!")
  }else if( !(tissueSiteDetail %in% c("All", tissueSiteDetailGTEx$tissueSiteDetail)) ){
    message("",paste0(c(paste0(1:nrow(qtl_tissue),". ",qtl_tissue$tissueSiteDetail)), collapse = "\n"))
    stop("Parameter \"tissueSiteDetail\" should be chosen from above tissue names!")
  }
  # convert tissueSiteDetail to tissueSiteDetailId:
  tissueSiteDetailId <- qtl_tissue[tissueSiteDetail, on ="tissueSiteDetail"]$tissueSiteDetailId

  message("== Querying eQTL associations from EBI API server:")
  eqtlAsso <- GTExdownload_assoAll(gene = gene, geneType = geneType,tissueSiteDetail= tissueSiteDetail, study= study )
  if( !exists("eqtlAsso") || is.null(eqtlAsso) || nrow(eqtlAsso)==0 ){
    stop("No eqtl associations were found for gene [", gene, "] in [", tissueSiteDetail, "] in [", study,"].")
  }else if( nrow(eqtlAsso)<=2 ){
    warning("Only ",nrow(eqtlAsso)," eqtl ", ifelse(nrow(eqtlAsso)==1,"association was","associations were")," detected, please extend the genome range using parameter \"coordinate\". " )
    return(NULL)
  }else{
    message("   Totally, ", nrow(eqtlAsso)," eQTL association was found in [", study, "] of gene [", gene,"] in [", tissueSiteDetail,"].")
    message("== Done")
  }

  # split coordinate from variantId:
  eqtlAsso <- cbind(eqtlAsso, data.table::rbindlist( lapply(eqtlAsso$variantId, function(x){ splitOut =stringr::str_split(x, stringr::fixed("_"))[[1]]; data.table(chrom=splitOut[1], pos=splitOut[2]) } )  ))
  eqtlAsso$logP <- ifelse(eqtlAsso$pValue==0,0, (-1*log(eqtlAsso$pValue, 10)))
  eqtlAsso$pos <- as.integer(eqtlAsso$pos)
  # order by logP desc:
  eqtlAsso <- eqtlAsso[order(-logP)]

  # title:
  plotTitle <- paste0(gene, " (", ifelse(stringr::str_detect(unique(eqtlAsso$chrom), stringr::regex("^chr")),unique(eqtlAsso$chrom), paste0("chr", unique(eqtlAsso$chrom))),
                      ":",
                      paste0(range(eqtlAsso$pos), collapse = "-")
                      ,")")

  # ylab and unit:
  posUnit <- "Bb"
  if( any(range(eqtlAsso$pos)>10^6)){
    eqtlAsso$pos <- eqtlAsso$pos/10^6
    posUnit <- "Mb"
  }else if( all(range(eqtlAsso$pos)<10^6) && all(range(eqtlAsso$pos)>10^3) ){
    eqtlAsso$pos <- eqtlAsso$pos/10^3
    posUnit <- "Kb"
  }else{
    posUnit <- "Bb"
  }
  yLab <-expression(-log["10"]("Pvalue"))

  # xlab:
  xLab <- paste0(ifelse(stringr::str_detect(unique(eqtlAsso$chrom), stringr::regex("^chr")),unique(eqtlAsso$chrom), paste0("chr", unique(eqtlAsso$chrom)))," (",posUnit,")")


  # color: 高亮最显著的那一个，并标注。
  eqtlAsso$colorP <- "grey"
  eqtlAsso$sizeP <- 1
  eqtlAsso[which(eqtlAsso$logP>=quantile(eqtlAsso$logP, seq(0,1,0.1))[10] & eqtlAsso$logP<quantile(eqtlAsso$logP, seq(0,1,0.1))[11] ),]$colorP <- "green"
  eqtlAsso[which(eqtlAsso$logP>=quantile(eqtlAsso$logP, seq(0,1,0.1))[10] & eqtlAsso$logP<quantile(eqtlAsso$logP, seq(0,1,0.1))[11] ),]$sizeP <- 1.01


  if( highlightSnp!="" && length(highlightSnp)==1 && !is.na(highlightSnp) && highlightSnp %in% eqtlAsso$snpId){
    eqtlAsso[snpId==highlightSnp,]$colorP <- "red"
    eqtlAsso[snpId==highlightSnp,]$sizeP <-  1.1
  }else{
    eqtlAsso[1,]$colorP <- "red"
    eqtlAsso[1,]$sizeP <- 1.1
  }

  if( requireNamespace("ggplot2") ){
    p <- ggplot(eqtlAsso)+
      geom_point(aes(x=pos, y=logP, color=colorP, size=sizeP))+
      scale_color_manual(breaks = eqtlAsso$colorP, values = eqtlAsso$colorP)+
      # scale_size_manual(breaks = as.character(eqtlAsso$sizeP), values =  as.numeric(eqtlAsso$sizeP) )+
      # geom_text(aes(x=pos, y=logP, label=snpId ))+
      geom_label_repel(data=eqtlAsso[1:1,], aes(x=pos, y=logP, label=snpId) )+
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
  }

  p


}






