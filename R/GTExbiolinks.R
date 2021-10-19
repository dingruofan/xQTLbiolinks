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

utils::globalVariables("tissueSiteDetailGTExv8")
utils::globalVariables("tissueSiteDetailGTExv7")
