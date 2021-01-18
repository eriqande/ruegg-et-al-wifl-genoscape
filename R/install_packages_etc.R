# just R code to install the packages needed for running the notebooks

# get the packages needed from CRAN
install.packages(
  c(
    "ade4",
    "car",
    "dggridR",
    "downloader",
    "ebirdst",
    "ecospat",
    "emdist",
    "fields",
    "geosphere",
    "ggspatial",
    "gridExtra",
    "hms",
    "igraph",
    "inlabru",
    "knitr",
    "lubridate",
    "mapplots",
    "MASS",
    "move",
    "raster",
    "RColorBrewer",
    "reticulate",
    "rgeos",
    "rmarkdown",
    "rnaturalearth",
    "rnaturalearthdata",
    "rubias",
    "rworldmap",
    "sessioninfo",
    "sf",
    "sp",
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
remotes::install_github("eriqande/TESS3_encho_sen")  # for special version of tess3r
