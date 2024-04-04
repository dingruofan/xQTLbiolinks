### xQTLbiolinks: a comprehensive and scalable tool for integrative analysis of molecular QTLs

**`xQTLbiolinks`** is a end-to-end bioinformatic tool for efficient mining and analyzing public and user-customized xQTLs data for the discovery of disease susceptibility genes. xQTLbiolinks consists of tailored functions that can be grouped into four modules: **Data retrieval**, **Pre-processing**, **Analysis** and **Visualization**.


Instructions, documentation, and tutorials can be found at [**here**](https://dingruofan.github.io/xQTLbiolinks/index.html).

<img src="https://private-user-images.githubusercontent.com/27127316/319440622-aa582137-7f9d-421a-b781-b3ec11f37d56.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MTIyMDQyMzQsIm5iZiI6MTcxMjIwMzkzNCwicGF0aCI6Ii8yNzEyNzMxNi8zMTk0NDA2MjItYWE1ODIxMzctN2Y5ZC00MjFhLWI3ODEtYjNlYzExZjM3ZDU2LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNDA0MDQlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjQwNDA0VDA0MTIxNFomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTAyMzk5NTk5YjNiMTc0NjI4NGNjZTk5Yjg5NjQ5MjA5YzIyZTE1YTZlNjZjMTgzZjkzNmZhNjIzMTJhMzU5YjEmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.1Dd-6Id95znyR-hCMPFZcF3vS5xrPoeEHQUcfXuHOfE" alt="Overview" width="100%" height="100%"/>

### Quick Start

1.  `xQTLbiolinks` can be installed and used on any operator systems supporting R. The latest version (v1.6.3) is also available at [GitHub repository](https://github.com/dingruofan/xQTLbiolinks/) and it can be installed through `devtools::install_github("dingruofan/xQTLbiolinks”)`. For more details, please refer to the instructions at **Installation** section below.
2. Find the [**Query and download**](https://dingruofan.github.io/xQTLbiolinks/articles/query_download.html) for xQTLs, gene, variant, tissue, sample and expressions.
3.  Find the [**Quick Start**](https://dingruofan.github.io/xQTLbiolinks/articles/Quick_start.html) for a quick application of colocalization analysis with xQTLbiolinks. Go through a whole [**Case study**](https://dingruofan.github.io/xQTLbiolinks/articles/Colocalization_analysis_with_xQTLbiolinks.html) of detection of casual vairants and genes in prostate cancer using `xQTLbiolinks`.
4.  The details and instructions of all functions implemented in xQTLbiolinks can be found [**here**](https://dingruofan.github.io/xQTLbiolinks/reference/index.html). Find more instructions with examples for visualizations [**here**](https://dingruofan.github.io/xQTLbiolinks/articles/visualization.html).

### Citation

If you find the xQTLbiolinks package or any of the source code in this repository useful for your work, please cite:

Ruofan Ding, Xudong Zou, Yangmei Qin, Lihai Gong, Hui Chen, Xuelian Ma, Shouhong Guang, Chen Yu, Gao Wang, Lei Li, **xQTLbiolinks: a comprehensive and scalable tool for integrative analysis of molecular QTLs**, *Briefings in Bioinformatics*, Volume 25, Issue 1, January 2024, bbad440,

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

------------------------------------------------------------------------
### Dependencies
#### The xQTLbiolinks package has the following dependencies: 
**R packages**: BiocGenerics, cowplot (>= 1.1.1), curl (>= 4.3.2), data.table (>= 1.14.2), DBI, SummarizedExperiment, GenomeInfoDb, GenomicFeatures, GenomicRanges, ggplot2 (>= 3.3.6), ggrepel, IRanges, jsonlite (>= 1.7.2), viridis, RMySQL, stringr (>= 1.4.0), utils (>= 4.0.3),VariantAnnotation, TxDb.Hsapiens.UCSC.hg38.knownGene, PupillometryR, coloc, hyprcoloc, knitr, rtracklayer, usethis, ggridges, CMplot, R.utils, ggforestplot.

### Installation
To install this R package, you will need to have required package `SummarizedExperiment` installed from Bioconductor with following command:
```
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("SummarizedExperiment") # For windows or linux
BiocManager::install("SummarizedExperiment",type="source") # For MAC
```

Once you have installed the required package, you can then install `xQTLbiolinks` from CRAN or github(recommended) using following command:
``` r
# Install from github to get the latest version.
if(!require("devtools")){install.packages("devtools")}
devtools::install_github("dingruofan/xQTLbiolinks")
```


