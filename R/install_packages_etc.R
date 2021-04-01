# just R code to install the packages needed for running the notebooks

# Note: on 2021-03-29 I found that dggrid had been removed from CRAN.
# But old versions (we installed 2.0.4) are archived and can be installed instead.

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
# These include the specific commit (ref) used
remotes::install_github("eriqande/genoscapeRtools", ref = "3e1dbfc0")
remotes::install_github("eriqande/snps2assays", ref = "8666f066")
remotes::install_github("eriqande/whoa", ref = "dddbeead")
remotes::install_github("eriqande/TESS3_encho_sen", ref = "2a0ce6c6")  # for special version of tess3r
