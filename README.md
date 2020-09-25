Ruegg et. al., WIFL Genoscape, Data and Code Repository
================
Kristen Ruegg, Eric C. Anderson, Marius Somvielle, Rachael A. Bay,
Christen Bossu, Mary Whitfield, Eben H. Paxton, Thomas B. Smith

  - [Overview](#overview)
  - [Preliminaries and Dependencies](#preliminaries-and-dependencies)
      - [Unix Programs](#unix-programs)
      - [Genomic resources](#genomic-resources)
      - [Geospatial data](#geospatial-data)
      - [Java jars](#java-jars)
      - [R Packages](#r-packages)
  - [RMarkdown Documents](#rmarkdown-documents)
      - [001-rad-seq-data-summaries-and-pca.Rmd](#rad-seq-data-summaries-and-pca.rmd)
      - [002-select-snps-for-assays-from-rad.Rmd](#select-snps-for-assays-from-rad.rmd)
  - [Literature Cited](#literature-cited)

**Last Updated:** 2020-09-22

# Overview

This repository includes code, data, and some intermediate results to
reproduce the results in Ruegg et al. (“Seasonal niche breadth predicts
population declines in a long-distance migratory bird”). You can get the
whole thing by cloning or downloading the repo from
<https://github.com/eriqande/ruegg-et-al-wifl-genoscape>.

If you are viewing this as the README of a GitHub repository, note that
you can read it in a somewhat friendlier format (i.e., with a table of
contents, etc.) on GitHub pages at:
<https://eriqande.github.io/ruegg-et-al-wifl-genoscape/>

RAD sequening and the production of BAMs and VCFs was described and
documented in Ruegg et al. (2018). The work described in this repository
uses the products from that work

The various steps in our analyses are described in a series of RMarkdown
documents with a fair bit of explanation/description of the procedures.
The few analyses that used cluster computing resources are described and
documented, but results from those analyses are included in
`stored_results` so that it is not necessary to have access to a cluster
to reproduce the results herein.

Running some of the RMarkdown documents requires a proper shell. These
were initially developed and run on a Mac.

Rmarkdown documents and scripts in the 000 series were developed
primarily by Eric C. Anderson. Rmarkdown documents and scripts in the
100 series were developed primarily by Marius Somveille.

All the RMarkdown documents or scripts can be evaluated en masse by
sourcing the R script `render-numbered-Rmds-and-scripts.R`. On a fairly
old mac laptop the run times for each are as follows:

    ## # A tibble: 2 x 2
    ##   File                                    `Running time in HH:MM:SS`
    ##   <chr>                                   <chr>                     
    ## 1 001-rad-seq-data-summaries-and-pca.Rmd  00:00:44                  
    ## 2 002-select-snps-for-assays-from-rad.Rmd 00:00:34

# Preliminaries and Dependencies

Below is a description of the needed dependencies. A script to install
all of them on one of our test clusters is at
`000-00-prepare-dependencies.sh` (NOT YET INCLUDED), but it may likely
need some tweaking on your system.

## Unix Programs

The following programs must be installed and available in the PATH.
Versions used on Eric’s laptop appear in parentheses. (NOT COMPLETE YET)

  - `bgzip`:
  - `samtools`:
  - `tabix`:

## Genomic resources

Create a directory in the top level of the repository, download the
willow flycatcher then modify the sequence names to be consistent with
the nomenclature used in this work and index the genome (that
nomenclature is already in in the FASTA file, but just not as the
sequence name so, we need to only shift some words around on each
sequence name line).

``` sh
mkdir genome
cd genome/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Empidonax_traillii/latest_assembly_versions/GCA_003031625.1_ASM303162v1/GCA_003031625.1_ASM303162v1_genomic.fna.gz
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Empidonax_traillii/latest_assembly_versions/GCA_003031625.1_ASM303162v1/GCA_003031625.1_ASM303162v1_assembly_report.txt

# now, change the sequence names and make sure that the DNA bases
# are all listed in UPPERCASE.
gzcat GCA_003031625.1_ASM303162v1_genomic.fna.gz | awk '
  /^>/ {
    tmp = substr($1, 2)
    newname = $6
    sub(/,/, "", newname)
    sub(/:/, "|", newname)
    $1 = ">" newname;
    $6 = tmp ",";
    print
    next
  }
  {print toupper($0)}
' > wifl-genome.fna
samtools faidx wifl-genome.fna
```

## Geospatial data

Make a directory in the top level of the repository called `geo-spatial`
and then download some large files from Natural Earth Data to there:

  - Natural Earth II with Shaded Relief, Water, and Drainages raster
    <https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/raster/NE2_HR_LC_SR_W_DR.zip>.
    Put or symlink the resulting directory, `NE2_HR_LC_SR_W_DR` into
    `geo-spatial`.
  - 10m-cultural-vectors, Admin 1 – States, Provinces, Download boundary
    lines:
    <https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces_lines.zip>.
    Put or symlink the resulting directory,
    `ne_10m_admin_1_states_provinces_lines` into `geo-spatial`.
  - Finally, get the coastlines:
    <https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip>,
    and put the resulting folder, `ne_10m_coastline` into `geo-spatial`.

## Climate data

Make a directory in the top level of the repository called `climate-data` with sub-directories `Temperature` and `Precipitation`. 

Download average temperature and precipitation at 2.5 minutes resolution from WorldClim at 
<https://worldclim.org/data/worldclim21.html>, and put the data into the corresponding sub-directories.


## Java jars

The following Java-based programs must be downloaded, and the paths to
their associated Jar files must be listed appropriately in the file
`script/java-jar-paths.R`. (NO JAVA JARS NEEDED AT THIS TIME…)

## R Packages

Packages must be downloaded from CRAN, BioConductor, and GitHub. The
following code, which is in `R/install_packages_etc.R`, will download
and install the necessary packages:

``` r
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
```

# RMarkdown Documents

The following RMarkdown documents should be evaluated in order. The
script `render-numbered-Rmds-and-scripts.R` will do that when run in the
top level of this repository. Some RMarkdown documents rely on outputs
from previous ones. Some of these RMarkdown documents include calls to
Unix utilities, so might not run on non-Unix or non-Linux architectures.

Outputs (figures, tables, R data objects, etc) from each RMarkdown
document are written to the `outputs/XXX` directories. Intermediate
files written during evaluation of each RMarkdown document are written
to the `intermediates/XXX` directories (THOUGH I MIGHT NOT USE
intermediates AT ALL). To facilitate working between the cluster and a
desktop/laptop, some outputs are written to the `stored_results/XXX`
directories which are version controlled and included in this repo.

Thumbnails of the figures and tables generated by each RMarkdown
document appear below (NOT YET). Additionally, at the top of each
section is a link to the compiled (HTML-version) of the RMarkdown
document on GitHub Pages.

## 001-rad-seq-data-summaries-and-pca.Rmd

(*Compiled RMarkdown HTML document on GitHub Pages:*
[001-rad-seq-data-summaries-and-pca.html](https://eriqande.github.io/ruegg-et-al-wifl-genoscape/001-rad-seq-data-summaries-and-pca.html))

Summarizing the RAD seq data, filtering out obvious paralogs and
monomorphic sites, then producing a PCA plot and estimates of pairwise
\(F_\mathrm{ST}\).

## 002-select-snps-for-assays-from-rad.Rmd

(*Compiled RMarkdown HTML document on GitHub Pages:*
[002-select-snps-for-assays-from-rad.html](https://eriqande.github.io/ruegg-et-al-wifl-genoscape/002-select-snps-for-assays-from-rad.html))

All the steps we went through to rank SNPs for ability to resolve
different populations/subspecies, and then the workflow for designing
Fluidigm assays from them.

# Literature Cited

<div id="refs" class="references">

<div id="ref-ruegg2018ecological">

Ruegg, Kristen, Rachael A Bay, Eric C Anderson, James F Saracco, Ryan J
Harrigan, Mary Whitfield, Eben H Paxton, and Thomas B Smith. 2018.
“Ecological Genomics Predicts Climate Vulnerability in an Endangered
Southwestern Songbird.” *Ecology Letters* 21 (7): 1085–96.

</div>

</div>
