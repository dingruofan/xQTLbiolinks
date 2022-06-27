### xQTLbiolinks: An R package aims to query, download, visualize and perform colocalization analyses between xQTL data and GWAS signals

**`xQTLbiolinks`** is a well-developed R package that enables users-customized query, extraction, and visualization of **molecular QTLs** (eQTLs and sQTLs) and **gene expression** data from public resources (e.g., GTEx) through the application programming interface (API) of [GTEx](https://gtexportal.org/home/api-docs/index.html) and [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/api-docs/).

xQTLbiolinks consists of tailored functions that can be grouped into four modules: **Query**, **Download**, **Analyze** and **Visualization**.

Instructions, documentation, and tutorials can be found at [**here**](https://dingruofan.github.io/xQTLbiolinks/index.html).

<img src="https://raw.githubusercontent.com/dingruofan/xQTLbiolinks/master/img/Overview.png" alt="Overview" width=100% height=100% />


### Quick Start
1. `xQTLbiolinks` has been successfully installed on Mac OS X, Linux, and Windows, using the `remotes` package to install directly from GitHub with command `remotes::install_github("dingruofan/xQTLbiolinks")`. Check more details in section "Installation" below.
2. Find the [**Full document**](https://dingruofan.github.io/xQTLbiolinks/articles/Quick_start.html) for a quick application of colocalization analysis with xQTLbiolinks.
3. Go through a whole [**Case study**](https://dingruofan.github.io/xQTLbiolinks/articles/Colocalization_analysis_with_xQTLbiolinks.html) of detection of casual vairants and genes in prostate cancer using `xQTLbiolinks`.
4. Then walk through these vignettes to learn more about xQTLbiolinks: [**Function Instruction**](https://dingruofan.github.io/xQTLbiolinks/reference/index.html) and [**Visualization of expression and xQTL**](https://dingruofan.github.io/xQTLbiolinks/articles/visualization.html).

### Citation
If you find the xQTLbiolinks package or any of the source code in this repository useful for your work, please cite:

>Ruofan Ding, Xudong Zou, Yangmei Qin, Gao Wang, Lei Li. **xQTLbiolinks: an R package for integrative analysis of xQTL data.** (submitted)

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

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

