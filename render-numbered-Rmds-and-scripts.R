

files <- c(
  "001-rad-seq-data-summaries-and-pca.Rmd",
  "002-select-snps-for-assays-from-rad.Rmd"
)

for(f in files) {

  outfile <- stringr::str_replace(f, "Rmd$", "stdout")
  errfile <- stringr::str_replace(f, "Rmd$", "stderr")

  message("")
  message("Processing         : " , f)
  message("Standard output to : ", outfile)
  message("Error output to    : ", errfile)
  message("")
  # a bunch of rigamoral to get it to render in an entirely fresh, clean R session.
  # See: https://github.com/rstudio/rmarkdown/issues/1204
  path <- callr::r(
    function(...) rmarkdown::render(...),
    args = list(input = f,
                output_format = "html_document",
                envir = globalenv()
    ),
    stdout = outfile,
    stderr = errfile
  )
}


# after generating those html files, this just hardlinks to them from the docs
# directory.  This is just a little thing to make it easy to commit them into the
# docs directory for serving up on GitHub pages, while gitignoring them in the main directory.
#system("cd docs; rm [012]??*.html; for i in ../[012]??*.html; do echo $i; ln $i $(basename $i);  done")


# Also, here let's copy the final figures over:
# message("Creating directory final-figs and copying figures to it")
# dir.create("final-figs", showWarnings = FALSE)
#
# file.copy("outputs/100/RoSA_figure1_map_with_inset.pdf", "final-figs/fig-01.pdf", overwrite = TRUE)
# file.copy("outputs/002/allele-frequencies.pdf", "final-figs/fig-02.pdf", overwrite = TRUE)
# file.copy("hand-edited-images/ultra-heatmap-chopped.pdf", "final-figs/fig-03.pdf", overwrite = TRUE)
# file.copy("outputs/105/RoSA_figure4_multipanel_estuary_gonadsi_fatness.pdf", "final-figs/fig-04.pdf", overwrite = TRUE)

