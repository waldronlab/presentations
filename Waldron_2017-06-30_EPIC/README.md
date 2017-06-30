# EPIC2017 "Functional Meta’omic Analyses for Microbial Communities"

## Preparation

You will need to install several pieces of software, in the following recommended order.

1. R (https://cran.r-project.org/)
    + Note that you should be running the **current** release version of R (3.4.0).  If you are running an older version, please upgrade it, otherwise you will likely have problems installing and running Bioconductor.
2. RStudio (https://www.rstudio.com/)
3. Bioconductor (http://bioconductor.org/install for instructions; installed from the R console in R or RStudio)
4. After Bioconductor is installed, install the following additional packages from the R console (R command line):

```
library(BiocInstaller)
useDevel(TRUE)
biocLite(c("phyloseq", "DESeq2", "gamlss", "knitr", "curatedMetagenomicData"))
```

All R packages from www.cran.r-project.org, www.bioconductor.org, and www.github.com can be installed using `biocLite()`, and this is the recommended way of installing all packages once you have installed Bioconductor. 

## Development vs release version

Active development occurs in the *development* version of Bioconductor, and the development version of `curatedMetagenomicData` currently provides more than twice as many samples as the release version. To use the release version of Bioconductor, skip the line `useDevel(TRUE)` above, or revert again to release by typing:

```
library(BiocInstaller)
useDevel(FALSE)
```
