% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{xQTLvisual_genesExp}
\alias{xQTLvisual_genesExp}
\title{Density plot of expression profiles for multiple genes}
\usage{
xQTLvisual_genesExp(
  genes,
  geneType = "auto",
  tissueSiteDetail = "",
  axis_text_size = 1.3,
  axis_title_size = 1.3,
  title_size = 1.4,
  xlab_text = "",
  ylab_text = "",
  title_text = "",
  color_map = "C"
)
}
\arguments{
\item{genes}{(character string or a character vector) gene symbol or gencode id (versioned or unversioned are both supported).}

\item{geneType}{(character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".}

\item{tissueSiteDetail}{(character) details of tissues in GTEx can be listed using \code{tissueSiteDetailGTExv8}}

\item{axis_text_size}{(numberic) text size of the axis labels}

\item{axis_title_size}{(numberic) text size of the axis title}

\item{title_size}{(numberic) text size of the title of the plot}

\item{xlab_text}{(character) Lable for x-axis}

\item{ylab_text}{(character) Lable for x-axis}

\item{title_text}{(character) Title of the plot}

\item{color_map}{A character string indicating the color map option to use. Eight options are available: "A", "B", "C", "D", "E", "F", "G", and "H"}
}
\value{
A ggplot object.
}
\description{
Density plot of expression profiles for multiple genes
}
\examples{
\donttest{
genes <- c("FNDC8", "S100Z", "AQP6", "AMOT", "C3orf38", "FOXL1", "COX11",
           "FCN3", "DDX58", "CFI", "MS4A18", "NUDT13", "HOXA4", "VSX1")
xQTLvisual_genesExp(genes, tissueSiteDetail="Lung")

genes <- c("ENSG00000073598.5","ENSG00000171643.13","ENSG00000086159.12","ENSG00000126016.15",
           "ENSG00000179021.9","ENSG00000176678.5","ENSG00000166260.10","ENSG00000142748.12",
           "ENSG00000107201.9","ENSG00000205403.12","ENSG00000214782.7","ENSG00000166321.13",
           "ENSG00000197576.13","ENSG00000100987.14")
xQTLvisual_genesExp(genes, geneType="gencodeId", tissueSiteDetail="Liver",  axis_text_size=1.6,axis_title_size=1.6, title_size=1.8, xlab_text="xxxx", ylab_text="yyyy", title_text ="tttt", color_map="E")
}
}
