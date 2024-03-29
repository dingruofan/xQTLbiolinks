% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{xQTLdownload_eqtlAllAsso}
\alias{xQTLdownload_eqtlAllAsso}
\title{Download summary statistics data for eQTLs with a specified gene, variant, tissue or study}
\usage{
xQTLdownload_eqtlAllAsso(
  gene = "",
  geneType = "auto",
  variantName = "",
  variantType = "auto",
  tissueLabel = "",
  study = "",
  recordPerChunk = 1000,
  withB37VariantId = FALSE,
  data_source = "eQTL_catalogue",
  API_version = "v1"
)
}
\arguments{
\item{gene}{(character) gene symbol or gencode id (versioned or unversioned are both supported).}

\item{geneType}{(character) options: "auto","geneSymbol" or "gencodeId". Default: "auto".}

\item{variantName}{(character) name of variant, dbsnp ID and variant id is supported, eg. "rs138420351" and "chr17_7796745_C_T_b38".}

\item{variantType}{(character) options: "auto", "snpId" or "variantId". Default: "auto".}

\item{tissueLabel}{(character) all supported tissues can be listed using "ebi_study_tissues"}

\item{study}{(character) name of studies can be listed using "ebi_study_tissues". If the study is null, use all studies (Default).}

\item{recordPerChunk}{(integer) number of records fetched per request (default: 1000).}

\item{withB37VariantId}{a logical value indicating whether to return the genome location(GTEx v7) of variants. Default: FALSE.}

\item{data_source}{"eQTL_catalogue"(default) or "liLab"}

\item{API_version}{"v1"(default) or "v2"(not working)}
}
\value{
A data.table object.
}
\description{
source of all eQTL associations is EBI eQTL category.
}
\examples{
\donttest{
# Download all eQTL associations of MLH1-rs13315355 pair in all tissues from all studies:
eqtlAsso <- xQTLdownload_eqtlAllAsso(gene="MLH1", variantName = "rs13315355")

# Download eQTL associations of gene ATP11B in CD4+ T cell from all supported studies:
geneAsso <- xQTLdownload_eqtlAllAsso(gene="MMP7",tissueLabel = "CD4+ T cell")

# Download eQTL associations of gene ATP11B in Muscle - Skeletal from GTEx_V8:
geneAsso <- xQTLdownload_eqtlAllAsso("ATP11B", tissueLabel="Muscle - Skeletal")

# Download eQTL associations of SNP rs11568818 in Muscle - Skeletal from GTEx_V8:
varAsso <- xQTLdownload_eqtlAllAsso(variantName="chr11_102530930_T_C_b38",
                                    tissueLabel="Muscle - Skeletal")

# Download all eQTL associations of SNP rs11568818 in all tissues from all supported studies.
varAsso <- xQTLdownload_eqtlAllAsso(variantName="rs11568818")
}
}
