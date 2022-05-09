# xQTLbiolinks: a R package aims to query, download, visual and integrate integrative analysis of GTEx data

By retrieving GTEx public-access data programmatically using the application programming interface (API) of GTEx and eQTL Catalogue, the functions provided in this package enable users to access molecular QTLs (eQTLs and sQTLs) and gene expressions data filtered by tissue, gene, variant or dataset, which assist in fast querying variants’ effect on gene expression in the process of disease studies; to uncover tissue-specific expressed genes; and to detect disease associated genes by conduct colocalization analyses of GWAS and eQTL signals; to perform data visualization (e.g. correlation plots of expression, boxplots of eQTL expression; and locuszom plots). xQTLbiolinks consists of functions that can be grouped into three main levels: Query, Download, Analysis and Visualization.

### Installation from GitHub

``` r
if(!require("remotes")){install.packages("remote")}
remotes::install_github("stronghoney/GTExbiolinks", 
                         host="api.gitee.com", 
                         auth_token="713472a1e408b426bc8da79acd5f7b12")
remotes::install_github("dingruofan/GTExbiolinks", 
                         auth_token="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")
```
