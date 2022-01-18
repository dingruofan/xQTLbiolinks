#' @title GTExbiolinks package
#' @description
#' The functions you're likely to need from \pkg{GTExbiolinks} is
#' Otherwise refer to the vignettes to see
#' how to format the documentation.
#'
#' @docType package
#' @name GTExbiolinks
NULL

#' @title  Tissue name and tissue id mapping of GTEx V8.
#' @description
#'  A dataset containing the 54 tissues' name and corresponding ID.
#' @docType data
#' @keywords internal
#' @name tissueSiteDetailGTExv8
#' @format A data frame with 54 rows and 2 variables
#' \describe{
#'   \item{tissueSiteDetail}{character string, tissue name}
#'   \item{tissueSiteDetailId}{character string, tissue id removding special character}
#' }
#' @source \url{https://gtexportal.org/home/}
NULL


#' @title Tissue name and tissue id mapping of GTEx V7.
#' @description
#'  A dataset containing the 53 tissues' name and corresponding ID.
#' @docType data
#' @keywords internal
#' @name tissueSiteDetailGTExv7
#' @format A data frame with 53 rows and 2 variables:
#' \describe{
#'   \item{tissueSiteDetail}{character string, tissue name}
#'   \item{tissueSiteDetailId}{character string, tissue id removding special character}
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

utils::globalVariables("tissueSiteDetailGTExv8")
utils::globalVariables("tissueSiteDetailGTExv7")
utils::globalVariables("sampleNum")
utils::globalVariables("gencodeGeneInfoAllGranges")





