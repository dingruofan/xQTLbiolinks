### xQTLbiolinks: a comprehensive and scalable tool for integrative analysis of molecular QTLs

**`xQTLbiolinks`** is a end-to-end bioinformatic tool for efficient mining and analyzing public and user-customized xQTLs data for the discovery of disease susceptibility genes. xQTLbiolinks consists of tailored functions that can be grouped into four modules: **Data retrieval**, **Pre-processing**, **Analysis** and **Visualization**.

Instructions, documentation, and tutorials can be found at [**here**](https://dingruofan.github.io/xQTLbiolinks/index.html).

<img src="https://raw.githubusercontent.com/dingruofan/xQTLbiolinks/master/img/Overview.png" alt="Overview" width="100%" height="100%"/>

### Quick Start

1.  `xQTLbiolinks` can be installed and used on any operator systems supporting R. Once the R (version 4.0 or later) is available, using `install.packages("xQTLbiolinks")` to install the steady version (v1.4.2) of `xQTLbiolinks`. The latest version (v1.4.3) is also available at [GitHub repository](https://github.com/dingruofan/xQTLbiolinks/) and it can be installed through `remotes::install_github("dingruofan/xQTLbiolinks”)`. For more details, please refer to the instructions at **Installation** section below.
2. Find the [**Query and download**](https://dingruofan.github.io/xQTLbiolinks/articles/query_download.html) for xQTLs, gene, variant, tissue, sample and expressions.
3.  Find the [**Quick Start**](https://dingruofan.github.io/xQTLbiolinks/articles/Quick_start.html) for a quick application of colocalization analysis with xQTLbiolinks. Go through a whole [**Case study**](https://dingruofan.github.io/xQTLbiolinks/articles/Colocalization_analysis_with_xQTLbiolinks.html) of detection of casual vairants and genes in prostate cancer using `xQTLbiolinks`.
4.  The details and instructions of all functions implemented in xQTLbiolinks can be found [**here**] (https://dingruofan.github.io/xQTLbiolinks/reference/index.html). Find more instructions with examples for visualizations [**here**](https://dingruofan.github.io/xQTLbiolinks/articles/visualization.html).

### Citation

If you find the xQTLbiolinks package or any of the source code in this repository useful for your work, please cite:

> Ruofan Ding, Xudong Zou, Yangmei Qin, Hui Chen, Xuelian Ma, Gao Wang, Lei Li. **xQTLbiolinks: a comprehensive and scalable tool for integrative analysis of molecular QTLs.** (submitted)

Institute of Systems and Physical Biology, Shenzhen Bay Laboratory, Shenzhen 518055, China

------------------------------------------------------------------------

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
if(!require("remotes")){install.packages("remotes")}
remotes::install_github("dingruofan/xQTLbiolinks")

# This command should automatically install any missing dependencies that are available from CRAN
install.packages("xQTLbiolinks")
```


