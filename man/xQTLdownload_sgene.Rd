% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GTExdownload.R
\name{xQTLdownload_sgene}
\alias{xQTLdownload_sgene}
\title{Fetch sGenes (sQTL Genes) from the specified dataset.}
\usage{
xQTLdownload_sgene(
  gene = "",
  geneType = "geneSymbol",
  datasetId = "gtex_v8",
  tissueSiteDetail = "",
  recordPerChunk = 2000
)
}
\arguments{
\item{gene}{A charater string of gene symbol, gencode id (versioned). Can be null.}

\item{geneType}{A character string. "geneSymbol"(default) or"gencodeId". Can be null.}

\item{datasetId}{A character string. only support "gtex_v8". Default: "gtex_v8".}

\item{tissueSiteDetail}{A character string.
Tissue must be chosen from the following tissue names:
\tabular{rrrrr}{
\strong{tissue name} \tab \strong{GTEx V8} \cr
Adipose - Subcutaneous \tab √ \cr
Adipose - Visceral (Omentum) \tab √ \cr
Adrenal Gland \tab √ \cr
Artery - Aorta \tab √ \cr
Artery - Coronary \tab √ \cr
Artery - Tibial \tab √ \cr
Bladder \tab √ \cr
Brain - Amygdala \tab √ \cr
Brain - Anterior cingulate cortex (BA24) \tab √ \cr
Brain - Caudate (basal ganglia) \tab √ \cr
Brain - Cerebellar Hemisphere \tab √ \cr
Brain - Cerebellum \tab √ \cr
Brain - Cortex \tab √ \cr
Brain - Frontal Cortex (BA9) \tab √ \cr
Brain - Hippocampus \tab √ \cr
Brain - Hypothalamus \tab √ \cr
Brain - Nucleus accumbens (basal ganglia) \tab √ \cr
Brain - Putamen (basal ganglia) \tab √ \cr
Brain - Spinal cord (cervical c-1) \tab √ \cr
Brain - Substantia nigra \tab √ \cr
Breast - Mammary Tissue \tab √ \cr
Cells - Cultured fibroblasts \tab √ \cr
Cells - EBV-transformed lymphocytes \tab √ \cr
Cells - Transformed fibroblasts \tab x \cr
Cervix - Ectocervix \tab √ \cr
Cervix - Endocervix \tab √ \cr
Colon - Sigmoid \tab √ \cr
Colon - Transverse \tab √ \cr
Esophagus - Gastroesophageal Junction \tab √ \cr
Esophagus - Mucosa \tab √ \cr
Esophagus - Muscularis \tab √ \cr
Fallopian Tube \tab √ \cr
Heart - Atrial Appendage \tab √ \cr
Heart - Left Ventricle \tab √ \cr
Kidney - Cortex \tab √ \cr
Kidney - Medulla \tab √ \cr
Liver \tab √ \cr
Lung \tab √ \cr
Minor Salivary Gland \tab √ \cr
Muscle - Skeletal \tab √ \cr
Nerve - Tibial \tab √ \cr
Ovary \tab √ \cr
Pancreas \tab √ \cr
Pituitary \tab √ \cr
Prostate \tab √ \cr
Skin - Not Sun Exposed (Suprapubic) \tab √ \cr
Skin - Sun Exposed (Lower leg) \tab √ \cr
Small Intestine - Terminal Ileum \tab √ \cr
Spleen \tab √ \cr
Stomach \tab √ \cr
Testis \tab √ \cr
Thyroid \tab √ \cr
Uterus \tab √ \cr
Vagina \tab √ \cr
Whole Blood \tab √ \cr
}}

\item{recordPerChunk}{A integer value (1-2000). number of records fetched per request (default: 2000).}
}
\value{
a data.table
}
\description{
sGenes are genes that have at least one significant sQTL acting upon them. Results may be filtered by tissue.
}
\examples{
\donttest{
 sGeneInfoAlltissue <- xQTLdownload_sgene()
 sGeneInfo <- xQTLdownload_sgene("TP53")
 sGeneInfo <- xQTLdownload_sgene(tissueSiteDetail="Lung", recordPerChunk=2000)
 sGeneInfo <- xQTLdownload_sgene("ENSG00000141510.16", geneType="gencodeId", tissueSiteDetail="Lung")
 sGeneInfo <- xQTLdownload_sgene("DDX11", tissueSiteDetail="Artery - Tibial" )
}
}