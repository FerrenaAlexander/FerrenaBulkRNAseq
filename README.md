
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FerrenaBulkRNAseq

<!-- badges: start -->
<!-- badges: end -->

This package goes through bulkRNAseq data anlysis with helpful basic
functions for comparative Differential Expression (DE) analysis and Gene
Set Enrichment Analysis(GSEA). DE testing uses DESeq2 and GSEA uses
FGSEA. Plotting is a big part of this package and relies on ggplot2 and
the general tidyverse suite.

## Installation

You can install from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("FerrenaAlexander/FerrenaBulkRNAseq")
```

The dependency package `fgsea` is a bioconductor package which seems to
require compilation to install, so if prompted you should install that
with compilation.

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(FerrenaBulkRNAseq) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
