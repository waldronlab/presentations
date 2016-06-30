# EPIC2016 "Functional Meta’omic Analyses for Microbial Communities"

## Preparation

You will need to install several pieces of software, in the following recommended order.

1. R (https://cran.r-project.org/)
    + Note that you should be running the **current** release version of R.  If you are running an older version, please upgrade it, otherwise you will likely have problems installing and running Bioconductor.
2. RStudio (https://www.rstudio.com/)
3. Bioconductor (http://bioconductor.org/install for instructions; installed from the R console in R or RStudio)
4. After Bioconductor is installed, install the following additional packages from the R console (R command line):

```
library(BiocInstaller)
biocLite(c("phyloseq", "DESeq2", "gamlss", "knitr"))
```

Note that all packages from both the CRAN and Bioconductor libraries can be installed in this way using `biocLite()`.  This is the recommended way of installing all packages once you have installed Bioconductor.

## Files

* `data/` datasets output from MetaPhlAn2:
    - `763577454.tsv`: HMP dataset
    - `Candela_Africa_metadat.txt`: participant metadata for Candela dataset
    - `Candela_Africa_stool.txt`: OTU data for Candela dataset
* `Waldron_stats.Rmd`: "R markdown" source file
* `Waldron_stats.R`: R code only

## Presentation

A built version of these slides is available at: http://rpubs.com/lwaldron/EPIC2016_Waldron