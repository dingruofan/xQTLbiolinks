## xQTLbiolinks: a R package aims to query, download, visual and perform integrative analyses of xQTL data

- **`xQTLbiolinks`** is an R package that enable users to access **molecular QTLs** (eQTLs and sQTLs) and **gene expressions** data filtered by tissue, gene, variant or dataset, by retrieving GTEx public-access data programmatically using the application programming interface (API) of [GTEx](https://gtexportal.org/home/api-docs) and [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/api-docs/).
- xQTLbiolinks consists of functions that can be grouped into four main levels: **Query**, **Download**, **Analyze** and **Visualize**.
<img src="img/b1ba7eb80950d93d626cb12acf4c54f5.png" alt="Overview" width=100% height=100% />
<br/>

### Quick Start

1. Install xQTLbiolinks from Github: `remotes::install_github("dingruofan/xQTLbiolinks")`. For more detailed installation instructions, see below.
2. See the [**Manual**](https://github.com/dingruofan/xQTLbiolinks/wiki/Colocalization-analysis-with-xQTLbiolinks) for a quick application of colocalization analysis with xQTLbiolinks .
3. Then walk through these vignettes to learn more about xQTLbiolinks: [Function Instruction](https://github.com/dingruofan/xQTLbiolinks/wiki/Function-Instruction) and [visualization of expression and xQTL](https://github.com/dingruofan/xQTLbiolinks/wiki/Visualization-of-expression-and-xQTL).

### Citation
Ruofan Ding, Xudong Zou, Gao Wang, Lei Li. xQTLbiolinks: an R/Bioconductor package for integrative analysis of xQTL data. (submitted)

>Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

***

### Installation

Package `SummarizedExperiment` is required but reposited in Biocouductor, instead of CRAN. So we first install `SummarizedExperiment` from Biocouductor, then install xQTLbiolinks from Github.

```r
# install required bioconductor package SummarizedExperiment:
if (!require("BiocManager")){install.packages("BiocManager")}
BiocManager::install("SummarizedExperiment")

# This command should automatically install any missing dependencies that are available from CRAN and GitHub.
if(!require("remotes")){install.packages("remotes")}
remotes::install_github("dingruofan/xQTLbiolinks")
```

