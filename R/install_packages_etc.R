# just R code to install the packages needed for running the notebooks

# get the packages needed from CRAN
install.packages(
  c(
    "hms",
    "knitr",
    "lubridate",
    "rmarkdown",
    "rubias",
    "sessioninfo",
    "sf",
    "tidyverse"
  ),
  repos = "http://cran.rstudio.com"
)

# Install packages from BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")



# Install Eric's Packages from GitHub
# UPDATE UPON FINAL RELEASE WITH THE ACTUAL COMMIT USED
remotes::install_github("eriqande/genoscapeRtools")
remotes::install_github("eriqande/snps2assays")
remotes::install_github("eriqande/whoa")
