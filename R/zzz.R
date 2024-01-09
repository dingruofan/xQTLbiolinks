##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://dingruofan.github.io/xQTLbiolinks/", "\n\n")

  # if (capabilities("libcurl")) {
  #   dl.method <- "libcurl"
  # } else {
  #   dl.method <- getOption("download.file.method", default = "auto")
  # }

  # options(clusterProfiler.download.method = dl.method)
  # options(timeout = max(300, getOption("timeout"))) # see ?download.file

  citation <- paste0("If you use ", pkgname, " in published research, please cite:\n\n",
                     "Ruofan Ding, Xudong Zou, Yangmei Qin, Lihai Gong, Hui Chen, Xuelian Ma, Shouhong Guang, Chen Yu, Gao Wang, Lei Li, xQTLbiolinks: a comprehensive and scalable tool for integrative analysis of molecular QTLs, Briefings in Bioinformatics, Volume 25, Issue 1, January 2024, bbad440",
                     "")

  packageStartupMessage(paste0(msg, citation))
}
