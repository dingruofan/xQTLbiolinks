#' @title Plot normalized expression among genotypes.
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
#'  expEqtl <- GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53",
#'                                tissueSiteDetail="Esophagus - Mucosa")
#'  expEqtl <- GTExvisual_eqtlExp(variantName="rs78378222", gene ="TP53",
#'                                tissueSiteDetail="Lung")
#'  expEqtl <- GTExvisual_eqtlExp(variantName="rs3778754", gene ="IRF5",
#'                                tissueSiteDetail="Whole Blood")
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
#' @param highlightSnp SNP of interst
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'. Default: 'ALL'.
#' @param recordPerChunk Number of records per fetch.
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
#' gwasFile <- tempfile(pattern = "file")
#' gwasURL <- "http://bioinfo.szbl.ac.cn/finalColoc/tmp/gwasFile/gwasChr6Sub1.txt"
#' utils::download.file(gwasURL, destfile=gwasFile)
#' gwasDF <- data.table::fread(gwasFile, sep="\t", header=TRUE)
#' gwasDF <- gwasDF[,c("rsid","P")]
#' gwasPlot <- GTExvisual_eqtlGWAS(gwasDF, traitGene="RP11-385F7.1", tissueSiteDetail="Lung")
#'
#' }
GTExvisual_eqtlGWAS <- function(gwasDF, traitGene="", traitGeneType="geneSymbol", tissueSiteDetail="", study_id="gtex_v8", highlightSnp="", population="EUR",  recordPerChunk=300){
  .<-NULL
  pValue <-logP.gwas <-logP.eqtl <- variant <- logP <- pos <- colorP <- sizeP <- snpId <- NULL
  pointShape<- r2Cut <- SNP_B <- R2 <- distance_0 <- snpId_2 <- snpId_1 <- pValue.gwas <- pValue.etql <-NULL
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
  population1000G <- c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')
  # rename:
  names(gwasDF) <- c("snpId", "pValue")
  gwasDF <- gwasDF[order(snpId, pValue)][!duplicated(snpId),]
  data.table::setindex(gwasDF, snpId)

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
  eqtlAsso <- GTExdownload_assoAll(gene=geneInfo$gencodeId, geneType = "gencodeId", tissueSiteDetail = tissueSiteDetail, recordPerChunk=recordPerChunk )
  eqtlAsso <- cbind(eqtlAsso, data.table::rbindlist(lapply(eqtlAsso$variantId, function(x){ x=stringr::str_split(x, stringr::fixed("_"))[[1]];return(data.table::data.table(chrom=x[1],pos=x[2])) })))
  eqtlAsso$pos <- as.integer(eqtlAsso$pos)
  # if no eqtl got:
  if( nrow(eqtlAsso)==0 || !exists("eqtlAsso")){
    stop("No eqtl associations were found for ",traitGeneType,": [",traitGene,"].")
  }

  # merge and obtain info：
  gwas_eqtl <- merge(gwasDF, eqtlAsso, by="snpId", sort=FALSE, suffixes = c(".gwas",".eqtl"))
  if( nrow(gwas_eqtl)==0 || !exists("gwas_eqtl") ){
    message("No intersection of GWAS and eQTL dataset!")
    return(data.table::data.table())
  }
  gwas_eqtl$logP.gwas <- ifelse(gwas_eqtl$pValue.gwas==0,0, (-1*log(gwas_eqtl$pValue.gwas, 10)))
  gwas_eqtl$logP.eqtl <- ifelse(gwas_eqtl$pValue.eqtl==0,0, (-1*log(gwas_eqtl$pValue.eqtl, 10)))
  gwas_eqtl[,"distance_0" :=.(sqrt(logP.gwas^2+logP.eqtl^2))]

  ######## Get LD:
  P_chrom <- stringr::str_replace_all( unique(gwas_eqtl$chrom), stringr::fixed("chr"),"")
  snpLD <- data.table::data.table()
  if( exists("highlightSnp") && highlightSnp!="" && length(highlightSnp)==1 && !is.na(highlightSnp) && highlightSnp %in% gwas_eqtl$snpId ){
    # Get LD info of all population:
    if( length(population)==1 && population=="ALL" ){
      for(i in 1:length(population1000G)){
        message("  Retrieving LD information of ",i,"/",length(population1000G)," : ",population1000G[i])
        ldTmp <- retrieveLD(P_chrom, highlightSnp, population1000G[i])
        if( !exists("ldTmp")||is.null(ldTmp)||nrow(ldTmp)==0 ){
          next()
        }else{
          snpLD <- rbind(snpLD, ldTmp)
        }
        rm(ldTmp)
      }
    }else if( all( population %in% population1000G) ){
      for(i in 1:length(population)){
        message("  Retrieving LD information of ",i,"/",length(population)," : ", population[i])
        ldTmp <- retrieveLD( P_chrom, highlightSnp, population[i])
        if( !exists("ldTmp")||is.null(ldTmp)||nrow(ldTmp)==0 ){
          next()
        }else{
          snpLD <- rbind(snpLD, ldTmp)
        }
        rm(ldTmp)
      }
    }else{
      message("Parameter \"population\" should be chosen from \"AFR\", \"AMR\", \"EAS\", \"EUR\", or \"SAS\". ")
      return(NULL)
    }
  }else{
    highlightSnp <- gwas_eqtl[which.max(distance_0)]$snpId
    # Get LD info of all population:
    if( length(population)==1 && population=="ALL" ){
      for(i in 1:length(population1000G)){
        message("  Retrieving LD information of ",i,"/",length(population1000G)," : ",population1000G[i])
        ldTmp <- retrieveLD( P_chrom, highlightSnp, population1000G[i])
        if( !exists("ldTmp")||is.null(ldTmp)||nrow(ldTmp)==0 ){
          next()
        }else{
          snpLD <- rbind(snpLD, ldTmp)
        }
        rm(ldTmp)
      }
    }else if( all( population %in% population1000G) ){
      for(i in 1:length(population)){
        message("  Retrieving LD information of ",i,"/",length(population)," : ",population[i])
        ldTmp <- retrieveLD( P_chrom, highlightSnp, population[i])
        if( !exists("ldTmp")||is.null(ldTmp)||nrow(ldTmp)==0 ){
          next()
        }else{
          snpLD <- rbind(snpLD, ldTmp)
        }
        rm(ldTmp)
      }
    }else{
      message("Parameter \"population\" should be chosen from \"AFR\", \"AMR\", \"EAS\", \"EUR\", or \"SAS\". ")
      return(NULL)
    }
  }

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
  gwas_eqtl <- merge(gwas_eqtl, snpLD[,.(snpId=SNP_B, r2Cut)], by ="snpId", all.x=TRUE, sort=FALSE)
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
  plotTitle <- paste0(traitGene, " (", paste0("chr",P_chrom),
                      ":",
                      paste0(range(gwas_eqtl$pos), collapse = "-")
                      ,")")

  # Xlab and ylab:
  xLab <- expression(-log["10"]("Pvalue (eQTL)"))
  yLab <- expression(-log["10"]("Pvalue (GWAS)"))

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
      theme(axis.text.x=element_text(size=rel(1.3)),
            axis.title.x=element_text(size=rel(1.3)),
            axis.title.y=element_text(size=rel(1.3)),
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
  return(gwas_eqtl)
}


#' @title eql trait plot
#'
#' @param gene gene name/gencode ID
#' @param geneType "geneSymbol" or "gencodeId"
#' @param tissueSiteDetail gene name
#' @param highlightSnp highlighted SNP
#' @param study "gtex_v8"
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'. Default: 'ALL'.
#' @param recordPerChunk Number of records per fetch.
#' @import data.table
#' @import stringr
#' @import ggplot2
#' @import ggrepel
#' @return A data.table
#' @export
#'
#' @examples
#' \donttest{
#'  geneAsso <- GTExvisual_eqtlTrait( gene="RP11-385F7.1", tissueSiteDetail = "Lung", population="EUR")
#' }
GTExvisual_eqtlTrait <- function(gene="", geneType="geneSymbol", highlightSnp="", tissueSiteDetail="", study ="gtex_v8", population="EUR", recordPerChunk=300){
  .<-NULL
  pointShape <- r2Cut <- SNP_B <- R2 <- logP <- pos <- colorP <- sizeP <- snpId<-NULL
  # gene = "RP11-385F7.1"
  # geneType = "geneSymbol"
  # highlightSnp=""
  # tissueSiteDetail = "Lung"
  # study = "gtex_v8"
  # population="EUR"
  # recordPerChunk=300

  population1000G <- c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')
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
  eqtlAsso <- GTExdownload_assoAll(gene = gene, geneType = geneType,tissueSiteDetail= tissueSiteDetail, study= study, recordPerChunk = recordPerChunk )
  eqtlAsso <- eqtlAsso[!is.na(snpId)]
  if( !exists("eqtlAsso") || is.null(eqtlAsso) || nrow(eqtlAsso)==0 ){
    stop("No eqtl associations were found for gene [", gene, "] in [", tissueSiteDetail, "] in [", study,"].")
  }else if( nrow(eqtlAsso)<=2 ){
    warning("Only ",nrow(eqtlAsso)," eqtl ", ifelse(nrow(eqtlAsso)==1,"association was","associations were")," detected, please extend the genome range using parameter \"coordinate\". " )
    return(NULL)
  }else{
    message("   Totally, [", nrow(eqtlAsso),"] eQTL associations were obtained of gene [", gene,"] in [", tissueSiteDetail,"] in [", study, "].")
    message("== Done")
  }

  # split coordinate from variantId:
  eqtlAsso <- cbind(eqtlAsso, data.table::rbindlist( lapply(eqtlAsso$variantId, function(x){ splitOut =stringr::str_split(x, stringr::fixed("_"))[[1]]; data.table(chrom=splitOut[1], pos=as.integer(splitOut[2])) } )  ))
  eqtlAsso$logP <- ifelse(eqtlAsso$pValue==0,0, (-1*log(eqtlAsso$pValue, 10)))
  P_chrom <- stringr::str_remove_all( eqtlAsso[1,]$chrom, stringr::fixed("chr"))

  # Get LD:
  snpLD <- data.table()
  if( exists("highlightSnp") && highlightSnp!="" && length(highlightSnp)==1 && !is.na(highlightSnp) && highlightSnp %in% eqtlAsso$snpId ){
    # Get LD info of all population:
    if( length(population)==1 && ( population %in% population1000G) ){
      snpLD <- retrieveLD( P_chrom, highlightSnp, population)
    }else{
      message("Parameter \"population\" should be chosen from \"AFR\", \"AMR\", \"EAS\", \"EUR\", or \"SAS\". ")
      return(NULL)
    }
  }else{
    highlightSnp <- eqtlAsso[which.max(logP)]$snpId
    # Get LD info of all population:
    if( length(population)==1 && ( population %in% population1000G) ){
      snpLD <- retrieveLD( P_chrom, highlightSnp, population)
    }else{
      message("Parameter \"population\" should be chosen from \"AFR\", \"AMR\", \"EAS\", \"EUR\", or \"SAS\". ")
      return(NULL)
    }
  }

  # Set LD SNP color:
  if(nrow(snpLD)>0){
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
  eqtlAsso <- merge(eqtlAsso, snpLD[,.(snpId=SNP_B, r2Cut)], by ="snpId", all.x=TRUE, sort=FALSE)
  eqtlAsso[is.na(r2Cut),"r2Cut"]<- "(0.0-0.2]"
  eqtlAsso[snpId==highlightSnp,"r2Cut"]<- "(0.8-1.0]"

  # set size:
  eqtlAsso$sizeP <- "small"
  eqtlAsso[logP<max(eqtlAsso$logP), "sizeP"] <- "small"
  eqtlAsso[logP>=2 & logP< max(eqtlAsso$logP), "sizeP"] <- "middle"
  eqtlAsso[logP == max(eqtlAsso$logP), "sizeP"] <- "large"
  eqtlAsso[snpId == highlightSnp, "sizeP"] <- "most"

  # set color:
  colorDT <- data.table( r2Cut = as.character(cut(c(0.2,0.4,0.6,0.8,1),breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('(0.0-0.2]','(0.2-0.4]','(0.4-0.6]','(0.6-0.8]','(0.8-1.0]'), include.lowest=TRUE)),
                         pointColor= c("#9C8B88", "#e09351", "#df7e66", "#b75347", "#A40340"),
                         pointFill = c("#9C8B88", "#e09351", "#df7e66", "#b75347", "#096CFD"),
                         pointSize = c(1,1,2,2,2.5))
  colorDT <- merge(colorDT, unique(eqtlAsso[,.(r2Cut)]), by="r2Cut",all.x=TRUE)[order(-r2Cut)]
  # set shape:
  eqtlAsso$pointShape <- "normal"
  eqtlAsso[snpId==highlightSnp,"pointShape"] <- "highlight"
  eqtlAsso$pointShape <- as.factor(eqtlAsso$pointShape)
  eqtlAsso <- eqtlAsso[order(r2Cut, logP)]
  eqtlAsso <- rbind( eqtlAsso[snpId!=highlightSnp ], eqtlAsso[snpId==highlightSnp, ] )

  # title:
  plotTitle <- paste0(gene, " (", ifelse(stringr::str_detect(P_chrom, stringr::regex("^chr")),P_chrom, paste0("chr", P_chrom)),
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
  yLab <- expression(-log["10"]("Pvalue"))

  # xlab:
  xLab <- paste0(ifelse(stringr::str_detect(P_chrom, stringr::regex("^chr")),P_chrom, paste0("chr", P_chrom))," (",posUnit,")")

  if( require("ggplot2") && require(ggrepel) ){
    p <- ggplot(eqtlAsso)+
      geom_point(aes(x=pos, y=logP, fill=r2Cut, color=r2Cut, size=pointShape, shape=pointShape))+
      scale_size_manual(breaks = c('normal', "highlight"), values =  c(2,3)  )+
      scale_shape_manual(breaks = c('normal', "highlight"), values =  c(16,23) )+
      scale_color_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointColor) +
      scale_fill_manual(expression("R"^2),breaks=colorDT$r2Cut, labels = colorDT$r2Cut, values = colorDT$pointFill) +
      # geom_text(aes(x=pos, y=logP, label=snpId ))+
      geom_label_repel(data=eqtlAsso[snpId==highlightSnp,], aes(x=pos, y=logP, label=snpId) )+
      labs(title = plotTitle )+
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
  }
  return(eqtlAsso[,-c("r2Cut", "pointShape")])
}



#' @title LocusZoom plot
#'
#' @param dataFrame A data.frame or a data.table object. Two columns are required: "snps", a character list, using an rsID or chromosome coordinate (e.g. "chr7:24966446"); P-value.
#' @param highlightSnp Default is the snp that with lowest p-value.
#' @param population
#' @param range
#' @param token
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#'  eqtlAsso <- GTExdownload_assoAll("KIF15", tissueSiteDetail = "Liver")
#'  DF <- eqtlAsso[,.(snpId, pValue)]
#' }
GTExvisual_locusZoom <- function( DF , highlightSnp="", population="EUR", range="", token="9246d2db7917"){




  if( require(LDlinkR) ){
    LDinfo <- LDmatrix(snps = gwasEqtlInfo$rsid,
                       pop = "EUR", r2d = "r2",
                       token = '9246d2db7917',
                       file =FALSE)
    LDmatrix(snps = "rs9381563", pop = "EUR", token = '9246d2db7917', file =FALSE)


    rownames(LDinfo) <-LDinfo$RS_number
    LDinfo$RS_number <- NULL
    gwasEqtlInfo <- gwasEqtlInfo[rsid %in% colnames(LDinfo)]
  }
}

#############
# library(LDlinkR)
# LDexpress(snps = c("rs345", "rs456"),
#           pop = c("YRI", "CEU"),
#           tissue = c("ADI_SUB", "ADI_VIS_OME"),
#           r2d = "r2",
#           r2d_threshold = "0.1",
#           p_threshold = "0.1",
#           win_size = "500000",
#           token = '9246d2db7917'
# )

# a <- LDproxy("rs456", "EUR", "r2", token = "9246d2db7917")
# a1 <- retrieve_LD("7",snp = "rs456", population = "EUR")
#
# LDpair(var1 = "rs3", var2 = "rs4", pop = "YRI", token = "9246d2db7917")
# LDpop(var1 = "rs3", var2 = "rs4", pop = "YRI", r2d = "r2",  token = "9246d2db7917")
#
# a <- LDproxy("rs2285691", "EUR", "r2", token = "9246d2db7917")
# a1 <- retrieve_LD("7",snp = "rs2285691", population = "EUR")
#
# SNPchip()
# list_gtex_tissues()
# list_chips()
#
# SNPclip(c("rs3", "rs4", "rs148890987"), "YRI", "0.1", "0.01", token = "9246d2db7917")
#
# system.time(b <- LDmatrix(snps = a$Coord[1:300],
#                           pop = "EUR", r2d = "r2",
#                           token = '9246d2db7917',
#                           file =FALSE))


