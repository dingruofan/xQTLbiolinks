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
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'  # EQTL associatons of TP53:
#'  eqtlInfo <- GTExdownload_eqtlSig(gene = "TP53", tissueSiteDetail = "Esophagus - Mucosa")
#'  eqtlExp <- GTExdownload_eqtlExp(variantName = eqtlInfo$snpId, gene = eqtlInfo$geneSymbol,
#'                                  tissueSiteDetail = "Esophagus - Mucosa")
#'  GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Esophagus - Mucosa")
#'  GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53", tissueSiteDetail="Lung")
#' }
GTExvisual_eqtlExp <- function(variantName="", gene="", variantType="snpId", geneType="geneSymbol", tissueSiteDetail="", datasetId="gtex_v8" ){
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

  # check network:
  pingOut <- apiAdmin_ping()
  if( !(!is.null(pingOut) && pingOut==200) ){
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

  eqtlInfo <- cbind(eqtlInfo, data.table::rbindlist(lapply(eqtlInfo$variantId, function(x){var_tmp = stringr::str_split(x,stringr::fixed("_"))[[1]]; data.table(chrom=var_tmp[1], pos=var_tmp[2], ref=var_tmp[3], alt=var_tmp[4])  })))
  # replace genotypes with ref and alt:
  genoLable <- data.table(genotypes=0:2, genoLabels = c( ifelse(nchar(eqtlInfo$ref)==1, paste0(rep(eqtlInfo$ref,2),collapse = ""), "Ref"),
                                            ifelse(nchar(eqtlInfo$ref)==1, paste0(c(eqtlInfo$ref, eqtlInfo$alt),collapse = ""), "Het"),
                                            ifelse(nchar(eqtlInfo$ref)==1, paste0(rep(eqtlInfo$alt,2),collapse = ""), "Hom")) )
  genoLable$genoLabels <- factor(genoLable$genoLabels, levels = genoLable$genoLabels)
  genoLable <- merge(eqtlExp, genoLable, by ="genotypes")
  # for Pie:
  genoLabelPie <- data.table(table(genoLable$genoLabels))
  names(genoLabelPie) <- c("genoLabels", "labelNum")
  genoLable <- merge(genoLable, genoLabelPie, by="genoLabels", sort = FALSE)


  p <- ggplot(genoLable)+
    geom_boxplot( aes(x= genoLabels, y=normExp, fill=genoLabels), alpha=0.6)+
    geom_bar( aes(x=Num, y=genoLabels, fill=genoLabels),width = 1, stat = "identity") + coord_polar("y", start=0)
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
  print(p)
  return(p)
}

# 用户设定范围：1. SNP 上下游，2. gene 上下游 3. 自定义范围。
# 根据范围获得GWAS的SNP。
# 根据GWAS SNP获得 eqtl 信息。

# GWAS数据：1. snpId, 2.chrom 3. pos

#' @title significance analysis of GTEX and eQTL
#'
#' @param gwasDF
#' @param queryTerm
#' @param queryType A character string. for variant, "snpId", "variantId". For gene, "geneSymbol" or "gencodeId".
#' @param flankUp A whole integer. KB. Default: 100
#' @param flankDown A whole integer. KB. Default: 100
#' @param coordinate
#' @param datasetId
#' @import data.table
#' @import stringr
#' @return
#' @export
#'
#' @examples
GTExanalyze_eqtlGWAS <- function(gwasDF, queryTerm="", queryType="snpId", flankUp=100, flankDown=100, coordinate="chr1:1-300000", datasetId="gtex_v8"){
  # Default:

  gwasDF <- data.table::fread("../GWAS_White-Blood-Cell-Traits_Tajuddin_2016.txt", sep="\t", header=TRUE)
  gwasDF <- gwasDF[trait=="White-Blood-Cell-Counts",][,.(rsid, chr, POS)]
  queryTerm="MED24"
  queryType="geneSymbol"
  tissueSiteDetail="Lung"
  datasetId="gtex_v7"
  flankUp=100
  flankDown=100


  cutNum <- 100
  # rename:
  names(gwasDF) <- c("snpId", "chrom", "pos")


  # parameter check: snpId
  message("Checking parameter!")
  if( !all(unlist(lapply(gwasDF$snpId, function(x){ stringr::str_detect(x,stringr::regex("^rs[0-9]{1,30}[0-9]$")) }))) ){
    message("The first column of \"gwasDF\" must be snp ID, which start with \"rs\"! ")
    return(data.table::data.table())
  }
  # parameter check: chorm
  if(datasetId=="gtex_v7"){
    if( !all( unlist( lapply( unique(gwasDF$chrom), function(x){ stringr::str_detect(x, stringr::regex("^1[0-9]$|^2[0-3]$|^[1-9]$|^[xXyY]$")) })) ) ){
      message("For dataset gtex_v7, The second column of \"gwasDF\" must be chromosome, which must be chosen from \"1-23, x, y\" ")
      return(data.table::data.table())
    }
  }else if(datasetId=="gtex_v8"){
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
  }else if( !queryType %in% c("snpId", "variantId", "geneSymbol", "gencodeId") ){
    stop("Parameter \"queryType\" must be chosen from \"snpId\", \"variantId\", \"geneSymbol\", \"gencodeId\"")
  }
  if( is.null(queryTerm) || is.na(queryTerm) || length(queryTerm)!=1 || queryTerm=="" ){
    stop("Parameter \"queryTerm\" can not be NULL or NA!")
  }

  #
  if(queryType == "geneSymbol"){
    geneInfo <- GTExquery_gene(queryTerm,gencodeVersion = "v19")
  }
  # snp in range:
  if(geneInfo$strand=="+"){
    geneInfo$rangeStart=geneInfo$start - as.integer( flankUp*1000)
    geneInfo$rangeEnd=geneInfo$start + as.integer( flankDown*1000 )
    gwasDF <- gwasDF[pos>=geneInfo$rangeStart & pos<=geneInfo$rangeEnd]
  }else if(geneInfo$strand=="-"){
    geneInfo$rangeStart=geneInfo$start - as.integer( flankDown*1000)
    geneInfo$rangeEnd=geneInfo$start + as.integer( flankUp*1000 )
    gwasDF <- gwasDF[pos>=geneInfo$rangeStart & pos<=geneInfo$rangeEnd]
  }
  # query snp:
  ##### query with GTExquery_varPos: START
  # query pos in batch per 100 terms.
  var_tmpCut <- cbind(gwasDF, data.table::data.table( ID=1:nrow(gwasDF), cutF = as.character(cut(1:nrow(gwasDF),breaks=seq(0,nrow(gwasDF)+cutNum,cutNum) )) ))
  var_tmpCut <- var_tmpCut[,.(posLis=list(pos)),by=c("chrom","cutF")]
  # out:
  var_tmpGot <- data.table::data.table()
  for( ii in 1:nrow(var_tmpCut)){
    message("   Got varints: ",nrow(var_tmpGot)+length(unlist(var_tmpCut[ii,]$posLis)),"/",nrow(gwasDF))
    var_tmpGot <- rbind(var_tmpGot, GTExquery_varPos(var_tmpCut[ii,]$chrom, pos=as.integer(unlist(var_tmpCut[ii,]$posLis)), datasetId = datasetId))
  }
  var_tmp <- merge(var_tmp[,.(variantId)], var_tmpGot, by="variantId", all.x=TRUE)




  gsawVariant <- data.table::data.table( snpId = gwasDF$rsid, chrom=gwasDF$snp_pos, )




}
