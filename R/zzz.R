##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://github.com/lilab-bioinfo/xQTLbiolinks/", "\n\n")

  # if (capabilities("libcurl")) {
  #   dl.method <- "libcurl"
  # } else {
  #   dl.method <- getOption("download.file.method", default = "auto")
  # }

  # options(clusterProfiler.download.method = dl.method)
  # options(timeout = max(300, getOption("timeout"))) # see ?download.file

  citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",
                     "Ruofan Ding, Xudong Zou, Yangmei Qin, Gao Wang, Lei Li. xQTLbiolinks: an R package for integrative analysis of xQTL data.",
                     "")


  packageStartupMessage(paste0(msg, citation))
}
