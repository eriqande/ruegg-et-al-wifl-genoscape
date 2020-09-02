# just R code to install the packages needed for running the notebooks

# get the packages needed from CRAN
install.packages(
  c(
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



# Get package genoscapeRtools at https://github.com/eriqande/genoscapeRtools
# This is a package of useful functions by Eric C. Anderson. We
# used commit `XXXXXX`.  Get it like this:
# remotes::install_github("eriqande/genoscapeRtools", ref = "6972defd")
# UPDATE UPON FINAL RELEASE WITH THE ACTUAL COMMIT USED
remotes::install_github("eriqande/genoscapeRtools")
