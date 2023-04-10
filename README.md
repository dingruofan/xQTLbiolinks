### xQTLbiolinks: xQTLbiolinks: a comprehensive and scalable tool for integrative analysis of molecular QTLs

**`xQTLbiolinks`** is a well-developed R package that enables users-customized data retrieval, pre-processing, analysis, and data visualization of **molecular QTLs** (eQTLs and sQTLs) and **gene expression** data from public resources (e.g., GTEx) through the application programming interface (API) of [GTEx](https://gtexportal.org/home/api-docs/index.html) and [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/api-docs/).

xQTLbiolinks consists of tailored functions that can be grouped into four modules: **Data retrieval**, **Pre-processing**, **Analysis** and **Visualization**.

Instructions, documentation, and tutorials can be found at [**here**](https://dingruofan.github.io/xQTLbiolinks/index.html).

<img src="https://raw.githubusercontent.com/dingruofan/xQTLbiolinks/master/img/Overview.png" alt="Overview" width="100%" height="100%"/>

### Quick Start

1.  `xQTLbiolinks` has been successfully installed on Mac OS X, Linux, and Windows, using the `install.package` function to install directly from CRAN with command `install.package("xQTLbiolinks")`. We recommend installing the latest version from GitHub (`remotes::install_github("dingruofan/xQTLbiolinks")`) to access more advanced software features. Check more details in section "Installation" below.
2.  Find the [**Full document**](https://dingruofan.github.io/xQTLbiolinks/articles/Quick_start.html) for a quick application of colocalization analysis with xQTLbiolinks.
3.  Go through a whole [**Case study**](https://dingruofan.github.io/xQTLbiolinks/articles/Colocalization_analysis_with_xQTLbiolinks.html) of detection of casual vairants and genes in prostate cancer using `xQTLbiolinks`.
4.  Then walk through these vignettes to learn more about xQTLbiolinks: [**Function Instruction**](https://dingruofan.github.io/xQTLbiolinks/reference/index.html) and [**Visualization of expression and xQTL**](https://dingruofan.github.io/xQTLbiolinks/articles/visualization.html).

### Citation

If you find the xQTLbiolinks package or any of the source code in this repository useful for your work, please cite:

> Ruofan Ding, Xudong Zou, Yangmei Qin, Xuelian Ma, Gao Wang, Lei Li. **xQTLbiolinks: a comprehensive and scalable tool for integrative analysis of molecular QTLs.** (submitted)

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

------------------------------------------------------------------------

### Installation

1. We strongly recommend to install in the environment of R, instead of Rstudio.

2. Package `SummarizedExperiment` is required but reposited in Biocouductor. If you haven't installed it, please install `SummarizedExperiment` from Biocouductor with following command.

``` r
# install required bioconductor package SummarizedExperiment:
if (!require("BiocManager")){install.packages("BiocManager")}
BiocManager::install("SummarizedExperiment")
```

3. Install `xQTLbiolinks` from CRAN or github.

``` r
# This command should automatically install any missing dependencies that are available from CRAN
install.packages("xQTLbiolinks")

# Or install from github to get the latest version.
if(!require("remotes")){install.packages("remotes")}
remotes::install_github("dingruofan/xQTLbiolinks")
```


