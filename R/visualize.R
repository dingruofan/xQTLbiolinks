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
#'   dot num+++
#'  # EQTL associatons of TP53:
#'  eqtlInfo <- GTExdownload_eqtlSig(gene = "TP53", tissueSiteDetail = "Esophagus - Mucosa")
#'  eqtlInfo <- GTExdownload_eqtlSig(gene = "IRF5", tissueSiteDetail = "Esophagus - Mucosa")
#'  eqtlExp <- GTExdownload_eqtlExp(variantName = eqtlInfo$snpId, gene = eqtlInfo$geneSymbol,
#'                                  tissueSiteDetail = "Esophagus - Mucosa")
#'  GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Esophagus - Mucosa")
#'  GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Lung")
#'
#'  GTExvisual_eqtlExp(variantName="rs4728150", gene ="IRF5", tissueSiteDetail="Esophagus - Mucosa")
#'
#'  GTExvisual_eqtlExp(variantName="rs3778754", gene ="IRF5", tissueSiteDetail="Whole Blood")
#' }
GTExvisual_eqtlExp <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8" ){
  genoLabels <- normExp <- labelNum <- p <- NULL
  # variantName="rs78378222"
  # gene="TP53"
  # tissueSiteDetail="Esophagus - Mucosa"
  # datasetId="gtex_v8"
  # variantType="snpId"
  # geneType="geneSymbol"

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
    message("== Done")
  }

  message("== Querying expression from API server:")
  suppressMessages(eqtlExp <- GTExdownload_eqtlExp(variantName = eqtlInfo$snpId, gene = eqtlInfo$geneSymbol, tissueSiteDetail = tissueSiteDetail))
  if( !exists("eqtlExp") || is.null(eqtlExp) || nrow(eqtlExp)==0 ){
    stop("No expression profiles were found for gene [", gene, "] in ", tissueSiteDetail, " in ", datasetId,".")
  }else{
    message("== Done")
  }

  eqtlInfo <- cbind(eqtlInfo, data.table::rbindlist(lapply(eqtlInfo$variantId, function(x){var_tmp = stringr::str_split(x,stringr::fixed("_"))[[1]]; data.table::data.table(chrom=var_tmp[1], pos=var_tmp[2], ref=var_tmp[3], alt=var_tmp[4])  })))
  # replace genotypes with ref and alt:
  genoLable <- data.table::data.table(genotypes=0:2, genoLabels = c( ifelse(nchar(eqtlInfo$ref)==1, paste0(rep(eqtlInfo$ref,2),collapse = ""), "Ref"),
                                            ifelse(nchar(eqtlInfo$ref)==1, paste0(c(eqtlInfo$ref, eqtlInfo$alt),collapse = ""), "Het"),
                                            ifelse(nchar(eqtlInfo$ref)==1, paste0(rep(eqtlInfo$alt,2),collapse = ""), "Hom")) )
  genoLable$genoLabels <- factor(genoLable$genoLabels, levels = genoLable$genoLabels)
  genoLable <- merge(eqtlExp, genoLable, by ="genotypes")
  # for Pie:
  genoLabelPie <- data.table::data.table(table(genoLable$genoLabels))
  names(genoLabelPie) <- c("genoLabels", "labelNum")
  genoLabelPie$legends <- paste0(genoLabelPie$genoLabels, "(",genoLabelPie$labelNum,")")

  stopifnot(require(ggplot2))
  stopifnot(require(ggrepel))
  p <- ggplot(genoLable)+
    geom_boxplot( aes(x= genoLabels, y=normExp, fill=genoLabels), alpha=0.6)+
    # scale_fill_manual(values=c("green", "red"))+
    # scale_fill_brewer(palette = "Dark2")+
    theme_bw()+
    xlab("Genotypes")+
    ylab("Normalized expression")+
    theme(axis.text.x=element_text(size=rel(1.3)),
          axis.text.y = element_text(size=rel(1.3)),
          axis.title = element_text(size=rel(1.3)),
          legend.position = "none",
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid",
                                           colour ="white")
    )
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
  return(p)
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
#'  gwasDF <- data.table::fread("../GWAS_Type-2-Diabetes_Wood_2016.txt.gz", sep="\t", header=TRUE)
#'  gwasDF <- gwasDF[,.(rsid, chr, snp_pos, pvalue, effect_allele, non_effect_allele)]
#'  GTExanalyze_eqtlGWAS(gwasDF, queryTerm="SFTPD", queryType="geneSymbol",
#'                       tissueSiteDetail="Artery - Tibial", datasetId="gtex_v7")
#'  GTExanalyze_eqtlGWAS(gwasDF, queryTerm="rs12778583", queryType="snpId",
#'                       tissueSiteDetail="Thyroid", datasetId="gtex_v7")
#'  GTExanalyze_eqtlGWAS(gwasDF, queryTerm="12:122404966-124404966", queryType="coordinate",
#'                       tissueSiteDetail="Whole Blood", datasetId="gtex_v7")
#' }
GTExanalyze_eqtlGWAS <- function(gwasDF, queryTerm="", queryType="snpId", eqtlTrait="", eqtlType="", tissueSiteDetail="", flankUp=100, flankDown=100, datasetId="gtex_v8"){
  pos <- pValue.gwas <- pValue.etql <-NULL
  # gwasDF <- data.table::fread("../GWAS_White-Blood-Cell-Traits_Tajuddin_2016.txt", sep="\t", header=TRUE)
  # gwasDF <- data.table::fread("../GWAS_Type-2-Diabetes_Wood_2016.txt.gz", sep="\t", header=TRUE)
  # gwasDF <- gwasDF[trait=="White-Blood-Cell-Counts",][,.(rsid, chr, POS)]
  # gwasDF <- gwasDF[,.(rsid, chr, snp_pos, pvalue, effect_allele, non_effect_allele)]

  # queryTerm="ABCB9"
  # queryType="geneSymbol"
  # eqtlTrait=""
  # eqtlType=""
  # tissueSiteDetail="Whole Blood"
  # datasetId="gtex_v7"
  # flankUp=100
  # flankDown=100

  # queryTerm="rs7953894"
  # queryType ="snpId"

  if( eqtlTrait=="" && eqtlType==""){
    eqtlTrait <- queryTerm
    eqtlType <- queryType
  }

  cutNum <- 100
  gencodeVersion <- "v26"
  # rename:
  gwasDF <- gwasDF[,1:4]
  names(gwasDF)[1:4] <- c("snpId", "chrom", "pos","pValue")

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
  }
  # parameter check: gwasDF : chromesome
  if( datasetId=="gtex_v7"){
    gencodeVersion <- "v19"
    genomeBuild="GRCh37/hg19"
    if( !all( unlist( lapply( unique(gwasDF$chrom), function(x){ stringr::str_detect(x, stringr::regex("^1[0-9]$|^2[0-3]$|^[1-9]$|^[xXyY]$")) })) ) ){
      message("For dataset gtex_v7, The second column of \"gwasDF\" must be chromosome, which must be chosen from \"1-23, x, y\" ")
      return(data.table::data.table())
    }
  }else if( datasetId=="gtex_v8"){
    gencodeVersion <- "v26"
    genomeBuild="GRCh38/hg38"
    if( !all( unlist( lapply( unique(gwasDF$chrom), function(x){ stringr::str_detect(x, stringr::regex("^chr1[0-9]$|^chr2[0-3]$|^chr[1-9]$|^chr[xXyY]$")) })) ) ){
      message("For dataset gtex_v8, the second column of \"gwasDF\" must be chromosome, which must be chosen from \"chr1-chr23, chrx, chry\" ")
      return(data.table::data.table())
    }
  }else{
    message("Parameter \"datasetId\" must be chosen from \"gtex_v7\" and \"gtex_v8\"  ")
    return(data.table::data.table())
  }
  # parameter check: pos
  if( any(!is.wholenumber(gwasDF$pos)) || any(gwasDF$pos<=0) ){
    stop("The third column of \"gwasDF\" must be an positive integer!")
  }
  # parameter check: queryType
  if( is.null(queryType) || is.na(queryType) || length(queryType)!=1 ){
    stop("Parameter \"queryType\" can not be NULL or NA!")
  }else if( !queryType %in% c("snpId", "variantId", "geneSymbol", "gencodeId", "coordinate") ){
    stop("Parameter \"queryType\" must be chosen from \"snpId\", \"variantId\", \"geneSymbol\", \"gencodeId\", and \"coordinate\" ")
  }
  if( is.null(queryTerm) || is.na(queryTerm) || length(queryTerm)!=1 || queryTerm=="" ){
    stop("Parameter \"queryTerm\" can not be NULL or NA!")
  }
  message("== Done")

  #
  if( queryType %in% c("geneSymbol", "gencodeId") ){
    message("== Querying gene info from API server:")
    # check network:
    bestFetchMethod <- apiAdmin_ping()
    if( !exists("bestFetchMethod") || is.null(bestFetchMethod) ){
      message("Note: API server is busy or your network has latency, please try again later.")
      return(NULL)
    }
    suppressMessages( queryInfo <- GTExquery_gene(queryTerm, geneType = eqtlType, gencodeVersion = gencodeVersion) )
    if( exists("queryInfo") && nrow(queryInfo)>0 ){
      message("== Done")
    }else{
      message("Can not get gene info for [", queryTerm,"] in ", datasetId,".")
    }
    # Obtain snps in  queried range:
    if(queryInfo$strand=="+"){
      queryInfo$rangeStart=queryInfo$start - as.integer( flankUp*1000)
      queryInfo$rangeEnd=queryInfo$start + as.integer( flankDown*1000 )
    }else if(queryInfo$strand=="-"){
      queryInfo$rangeStart=queryInfo$start - as.integer( flankDown*1000)
      queryInfo$rangeEnd=queryInfo$start + as.integer( flankUp*1000 )
    }
    eqtlAsso <- GTExdownload_eqtlAll(gene=eqtlTrait, geneType = eqtlType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId )
  }else if( queryType %in% c("snpId", "variantId") ){
    queryInfo <- GTExquery_varId(queryTerm, variantType = eqtlType, datasetId = datasetId)
    queryInfo$rangeStart=queryInfo$pos - as.integer( flankUp*1000)
    queryInfo$rangeEnd=queryInfo$pos + as.integer( flankDown*1000 )
    # get all snps in a range.
    rangedVariants <- dbsnpQueryRange(ifelse(datasetId=="gtex_v8", queryInfo$chromosome, paste0("chr",queryInfo$chromosome)),
                                      startPos = queryInfo$rangeStart, endPos = queryInfo$rangeEnd, genomeBuild = genomeBuild )
    if(exists("rangedVariants") && nrow(rangedVariants)>1){
      message("== Totally ",  nrow(rangedVariants), " variants detected in genome coordinate: ", paste0(queryInfo$chromosome,":",queryInfo$rangeStart,"-",queryInfo$rangeEnd),".")
    }else{
      stop("No variants detected in genome coordinate: ", paste0(queryInfo$chromosome,":",queryInfo$rangeStart,"-",queryInfo$rangeEnd),".")
    }

    eqtlAsso <- GTExdownload_eqtlAll(gene=eqtlTrait, variantType = eqtlType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId )
    if(nrow(eqtlAsso) == 0 || !exists(eqtlAsso)){
      message("No significant associations were found for ",eqtlType," [",eqtlTrait,"] in Whole Blood in gtex_v7")
      return(NULL)
    }
  }else if ( queryType == "coordinate" ){
    message("Querying coordinate!")
    queryInfo <- stringr::str_split(queryTerm, stringr::regex("[:-]"))[[1]]
    if( is.null(queryInfo) || is.na(queryInfo) || length(queryInfo)!=3 ){
      message("For dataset gtex_v7, the format of specified genomic coordinate should be: \"chromesome:start-end\", like: \"chr1:1-300000\" ")
      return(NULL)
    }
    # check the first splited string:
    if( datasetId=="gtex_v7"){
      if( !stringr::str_detect(queryInfo[1], stringr::regex("^1[0-9]$|^2[0-3]$|^[1-9]$|^[xXyY]$")) ){
        message("For dataset gtex_v7, the format of specified genomic coordinate should be: \"chromesome:start-end\", like: \"1:1-300000\" ")
        message("For dataset gtex_v8, the format of specified genomic coordinate should be: \"chromesome:start-end\", like: \"chr1:1-300000\" ")
        return(NULL)
      }
    }else if( datasetId=="gtex_v8"){
      if( !stringr::str_detect(queryInfo[1], stringr::regex("^chr1[0-9]$|^chr2[0-3]$|^chr[1-9]$|^chr[xXyY]$")) ){
        message("For dataset gtex_v7, the format of specified genomic coordinate should be: \"chromesome:start-end\", like: \"1:1-300000\" ")
        message("For dataset gtex_v8, the format of specified genomic coordinate should be: \"chromesome:start-end\", like: \"chr1:1-300000\" ")
        return(NULL)
      }
    }
    # check the last splited string:
    if( any(!is.wholenumber(as.integer(queryInfo[2:3]))) || any(as.integer(queryInfo[2:3])<=0) || diff(as.integer(queryInfo[2:3]))<=0  ){
      message("For dataset gtex_v7, the format of specified genomic coordinate should be: \"chromesome:start-end\", like: \"1:1-300000\" ")
      message("For dataset gtex_v8, the format of specified genomic coordinate should be: \"chromesome:start-end\", like: \"chr1:1-300000\" ")
      return(NULL)
    }

    queryInfo <- data.table(chrom=queryInfo[1], start=queryInfo[2], end=queryInfo[3])
    queryInfo$rangeStart = queryInfo$start
    queryInfo$rangeEnd = queryInfo$end
    eqtlAsso <- GTExdownload_eqtlAll(variantName=eqtlTrait, variantType = eqtlType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId )
  }else{
    stop("Please enter the right parameter of \"queryType\". ")
  }
  # if no eqtl got:
  if(nrow(eqtlAsso)==0 || !exists("eqtlAsso")){
    stop("No eqtl associations were found for ",queryType,": [",queryTerm,"].")
  }
  gwasDFsub <- gwasDF[pos>=queryInfo$rangeStart & pos<=queryInfo$rangeEnd]
  if(nrow(gwasDFsub)==0 || !exists("gwasDFsub")){
    stop("Queried range is too small. please extend queried range.")
  }

  # query snp:
  ##### query with GTExquery_varPos: START
  # query pos in batch per 100 terms.
  # var_tmpCut <- cbind(gwasDFsub, data.table::data.table( ID=1:nrow(gwasDFsub), cutF = as.character(cut(1:nrow(gwasDFsub),breaks=seq(0,nrow(gwasDFsub)+cutNum,cutNum) )) ))
  # var_tmpCut <- var_tmpCut[,.(posLis=list(pos)),by=c("chrom","cutF")]
  # if(nrow(var_tmpCut)>100){
  #   message("Too many variants in queried range, please reset range!")
  # }
  # # out:
  # var_tmpGot <- data.table::data.table()
  # for( ii in 1:nrow(var_tmpCut)){
  #   message("   Got varints: ",nrow(var_tmpGot)+length(unlist(var_tmpCut[ii,]$posLis)),"/",nrow(gwasDFsub))
  #   var_tmpGot <- rbind(var_tmpGot, GTExquery_varPos(var_tmpCut[ii,]$chrom, pos=as.integer(unlist(var_tmpCut[ii,]$posLis)), datasetId = datasetId))
  # }
  # var_tmp <- merge(var_tmp[,.(variantId)], var_tmpGot, by="variantId", all.x=TRUE)

  gwas_eqtl <- merge(gwasDFsub, eqtlAsso, by="snpId", sort=FALSE, suffixes = c(".gwas",".etql"))
  if( nrow(gwas_eqtl)==0 || !exists("gwas_eqtl") ){
    message("No intersection of GWAS and eQTL dataset!")
    return(data.table::data.table())
  }

  stopifnot(require(ggplot2))
  stopifnot(require(ggrepel))
  p<- ggplot(gwas_eqtl)+
    geom_point( aes(x= (-1) * log(pValue.gwas,10), y = (-1)* log(pValue.etql,10)) )+
    theme_bw()+
    theme(axis.title.x=element_text(size=rel(1.3)),
          axis.title.y=element_text(size=rel(1.3)),
          legend.title = element_text(face="bold",size=rel(1.1)),
          legend.text = element_text(size=rel(1.1)),
          legend.position = "none",
          )
  print(p)
  return(p)

}


#' @title
#'
#' @param gene
#' @param geneType
#' @param coordinate
#' @param tissueSiteDetail
#' @param datasetId
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'  GTExvisual_eqtlTrait("ENSG00000112137.12", "gencodeId",tissueSiteDetail = "Adipose - Subcutaneous", datasetId="gtex_v7")
#' }
GTExvisual_eqtlTrait <- function(gene="", geneType="geneSymbol", coordinate="chr1:1-3", tissueSiteDetail="", datasetId="gtex_v8"){

  # gene = "ENSG00000112137.12"
  # geneType = "gencodeId"
  # tissueSiteDetail = "Adipose - Subcutaneous"
  # datasetId = "gtex_v7"

  eqtlOfgene <- GTExdownload_eqtlAll( gene=gene, geneType = geneType, tissueSiteDetail = tissueSiteDetail, datasetId = datasetId)
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
  eqtlOfgene$sizeP <- 1
  eqtlOfgene[which(eqtlOfgene$logP>=quantile(eqtlOfgene$logP, seq(0,1,0.1))[10] & eqtlOfgene$logP<quantile(eqtlOfgene$logP, seq(0,1,0.1))[11] ),]$colorP <- "green"
  eqtlOfgene[which(eqtlOfgene$logP>=quantile(eqtlOfgene$logP, seq(0,1,0.1))[10] & eqtlOfgene$logP<quantile(eqtlOfgene$logP, seq(0,1,0.1))[11] ),]$sizeP <- 2
  eqtlOfgene[1,]$colorP <- "red"
  eqtlOfgene[1,]$sizeP <- 3

  stopifnot(require(ggplot2))
  stopifnot(require(ggrepel))
  p <- ggplot(eqtlOfgene)+
    geom_point(aes(x=pos, y=logP, color=colorP, size=colorP))+
    scale_color_manual(breaks = eqtlOfgene$colorP, values = eqtlOfgene$colorP)+
    scale_size_manual(breaks = eqtlOfgene$sizeP, values =  rel(as.numeric(eqtlOfgene$sizeP)) )+
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
    )
  print(p)
  p


1


}






