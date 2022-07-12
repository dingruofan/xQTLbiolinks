#' @title xQTLbiolinks package
#' @description
#' the functions provided in this \pkg{xQTLbiolinks} enable users to access molecular QTLs (eQTLs and sQTLs) and gene expressions data filtered by tissue, gene, variant or dataset.
#'
#' @docType package
#' @name xQTLbiolinks
NULL

#' @title  Tissue name and tissue id mapping of GTEx V8.
#' @description
#'  A dataset containing the 54 tissues' name and corresponding ID of GTEx V8.
#' @docType data
#' @keywords internal
#' @name tissueSiteDetailGTExv8
#' @format A data frame with 54 rows and 2 variables
#' \describe{
#'   \item{tissueSiteDetail}{character string, tissue name}
#'   \item{tissueSiteDetailId}{character string, tissue id removding special character}
#'   \item{tissueSite}{character string, whole tissue name}
#'   \item{colorHex}{character string, color of the tissue in GTEx}
#' }
#'  Compared with GTEx v7
#' \tabular{rrrrr}{
#'   \strong{Tissue name} \tab \strong{GTEx V8} \tab \strong{GTEx V7} \cr
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
#' @source \url{https://gtexportal.org/home/}
NULL


#' @title Tissue name and tissue id mapping of GTEx V7.
#' @description
#'  A dataset containing the 53 tissues' name and corresponding ID of GTEx V7.
#' @docType data
#' @keywords internal
#' @name tissueSiteDetailGTExv7
#' @format A data frame with 53 rows and 2 variables:
#' \describe{
#'   \item{tissueSiteDetail}{character string, tissue name}
#'   \item{tissueSiteDetailId}{character string, tissue id removding special character}
#'   \item{tissueSite}{character string, whole tissue name}
#'   \item{colorHex}{character string, color of the tissue in GTEx}
#' }
#'  Compared with GTEx v8
#' \tabular{rrrrr}{
#'   \strong{Tissue name} \tab \strong{GTEx V8} \tab \strong{GTEx V7} \cr
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
#' @source \url{https://gtexportal.org/home/}
NULL

#' @title samples used in GTEx eQTL analysis.
#' @description
#'  A dataset containing the 49 tissues' name and corresponding sample number.
#' @docType data
#' @keywords internal
#' @name sampleNum
#' @format A data frame with 49 rows and 2 variables:
#' \describe{
#'   \item{tissueSiteDetailId}{character string, tissue name}
#'   \item{sampleNum}{integer}
#' }
#' @source \url{https://gtexportal.org/home/}
NULL

#' @title Gene annotations (chr1-chr22).
#' @description
#'  A dataset containing the gene information.
#' @docType data
#' @keywords internal
#' @name gencodeGeneInfoAllGranges
#' @format A GRanges object
#' \describe{
#'   \item{seqnames}{character string, chromosome}
#'   \item{ranges}{Iranges, gene location of v26}
#'   \item{strand}{character string, strand}
#'   \item{rangesV19}{Iranges, gene location of v19}
#'   \item{gencodeId}{character, gencode id}
#' }
#' @source \url{https://gtexportal.org/home/}
NULL

#' @title Gene types in GTEx V8 and V7
#' @description
#'  A dataset containing the classification of genes
#' @docType data
#' @keywords internal
#' @name gencodeGenetype
#' @format A list
#' \describe{
#'   \item{V26}{A character vector}
#'   \item{V19}{A character vector}
#' }
#' @source \url{https://gtexportal.org/home/}
NULL

#' @title data for vignette
#' @description
#'  head 6 rows of gwasDF of GLGC_CG0276_result.txt(GCST006085)
#' @docType data
#' @keywords internal
#' @name example_Coloc_gwasDF
#' @format A Data.table
#' \describe{
#'   \item{rsid}{A character vector}
#'   \item{chr}{A character vector}
#'   \item{position}{A integer vector}
#'   \item{pValue}{A numeric vector}
#'   \item{maf}{A numeric vector}
#'   \item{beta}{A numeric vector}
#'   \item{se}{A numeric vector}
#' }
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006085/harmonised/29892016-GCST006085-EFO_0001663-build37.f.tsv.gz}
NULL

#' @title data for vignette
#' @description
#'  head 6 rows of sentinelSnpDF of xQTLanalyze_getSentinelSnp
#' @docType data
#' @keywords internal
#' @name example_Coloc_sentinelSNP
#' @format A Data.table
#' \describe{
#'   \item{rsid}{A character vector}
#'   \item{chr}{A character vector}
#'   \item{position}{A integer vector}
#'   \item{pValue}{A numeric vector}
#'   \item{maf}{A numeric vector}
#'   \item{beta}{A numeric vector}
#'   \item{se}{A numeric vector}
#' }
#' @source xQTLanalyze_getSentinelSnp
NULL

#' @title data for vignette
#' @description
#'  head 6 rows of traitsAll of xQTLanalyze_getTraits
#' @docType data
#' @keywords internal
#' @name example_Coloc_traitsAll
#' @format A Data.table
#' \describe{
#'   \item{chromosome}{A character vector}
#'   \item{geneStart}{A integer vector}
#'   \item{geneEnd}{A integer vector}
#'   \item{geneStrand}{A character vector}
#'   \item{geneSymbol}{A character vector}
#'   \item{gencodeId}{A character vector}
#'   \item{rsid}{A character vector}
#'   \item{position}{A integer vector}
#'   \item{pValue}{A numeric vector}
#'   \item{maf}{A numeric vector}
#' }
#' @source xQTLanalyze_getTraits
NULL

#' @title data for vignette
#' @description
#'  head 6 rows of traitsAll of xQTLanalyze_coloc
#' @docType data
#' @keywords internal
#' @name example_Coloc_colocResultAll
#' @format A Data.table
#' \describe{
#'   \item{traitGene}{A character vector}
#'   \item{nsnps}{A character vector}
#'   \item{PP.H0.abf}{A numeric vector}
#'   \item{PP.H1.abf}{A numeric vector}
#'   \item{PP.H2.abf}{A numeric vector}
#'   \item{PP.H3.abf}{A numeric vector}
#'   \item{PP.H4.abf}{A numeric vector}
#'   \item{candidate_snp}{A character vector}
#'   \item{hypr_posterior}{A numeric vector}
#'   \item{hypr_regional_prob}{A numeric vector}
#'   \item{hypr_candidate_snp}{A character vector}
#'   \item{hypr_posterior_explainedBySnp}{A numeric vector}
#' }
#' @source xQTLanalyze_coloc
NULL

#' @title data for vignette
#' @description
#'  head 5 rows of hyprcoloc results from xQTLanalyze_coloc
#' @docType data
#' @keywords internal
#' @name example_Coloc_hyprcolocResultAll
#' @format A Data.table
#' \describe{
#'   \item{traitGene}{A character vector}
#'   \item{posterior_prob}{A numeric vector}
#'   \item{regional_prob}{A numeric vector}
#'   \item{candidate_snp}{A character vector}
#'   \item{posterior_explained_by_snp}{A numeric vector}
#' }
#' @source xQTLanalyze_coloc
NULL

#' @title data for vignette
#' @description
#'  head 6 rows of traitsAll of significant results from xQTLanalyze_coloc
#' @docType data
#' @keywords internal
#' @name example_Coloc_colocResultsig
#' @format A Data.table
#' \describe{
#'   \item{traitGene}{A character vector}
#'   \item{nsnps}{A character vector}
#'   \item{PP.H0.abf}{A numeric vector}
#'   \item{PP.H1.abf}{A numeric vector}
#'   \item{PP.H2.abf}{A numeric vector}
#'   \item{PP.H3.abf}{A numeric vector}
#'   \item{PP.H4.abf}{A numeric vector}
#'   \item{candidate_snp}{A character vector}
#'   \item{SNP.PP.H4}{A numeric vector}
#'   \item{hypr_posterior}{A numeric vector}
#'   \item{hypr_regional_prob}{A numeric vector}
#'   \item{hypr_candidate_snp}{A character vector}
#'   \item{hypr_posterior_explainedBySnp}{A numeric vector}
#' }
#' @source xQTLanalyze_coloc
NULL

#' @title data for vignette
#' @description
#'  a dataset containing all study-tissues mapping relationships.
#' @docType data
#' @keywords internal
#' @name ebi_study_tissues
#' @format A Data.table
#' \describe{
#'   \item{study_accession}{A character vector}
#'   \item{tissue_label}{A character vector}
#'   \item{tissue}{A character vector}
#' }
#' @source \url{https://www.ebi.ac.uk/eqtl/api-docs/#accessing-the-api} or EBIquery_allTerm
NULL

utils::globalVariables("tissueSiteDetailGTExv8")
utils::globalVariables("tissueSiteDetailGTExv7")
utils::globalVariables("sampleNum")
utils::globalVariables("gencodeGeneInfoAllGranges")
utils::globalVariables("gencodeGenetype")
# utils::globalVariables("ebi_qtl_groups")
utils::globalVariables("ebi_study_tissues")

utils::globalVariables("example_Coloc_gwasDF")
utils::globalVariables("example_Coloc_sentinelSNP")
utils::globalVariables("example_Coloc_traitsAll")
utils::globalVariables("example_Coloc_colocResultAll")
utils::globalVariables("example_Coloc_colocResultsig")
utils::globalVariables("example_Coloc_hyprcolocResultAll")

