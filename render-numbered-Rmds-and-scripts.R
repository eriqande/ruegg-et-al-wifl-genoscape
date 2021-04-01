

files <- c(
  "001-rad-seq-data-summaries-and-pca.Rmd",
  "002-select-snps-for-assays-from-rad.Rmd",
  "003-process-5-plates-of-birds-at-192-fluidigm-snps.Rmd",
  "004-combine-RAD-and-fludigm-breeders-and-run-STRUCTURE.Rmd",
  "005-choosing-96-SNPs-from-amongst-the-179.Rmd",
  "006-make-the-genoscape.Rmd",
  "007-process-fluidigm-plates-for-wintering-birds.Rmd",
  "008-rubias-self-assign-and-wintering-birds.Rmd"
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

# and, finally, run 101
path <- callr::r(
  function(...) source(...),
  args = list(file = "101-niche_tracking_analysis.R"
  ),
  stdout = "101-niche_tracking_analysis.stdout",
  stderr = "101-niche_tracking_analysis.stderr"
)



# after generating those html files, this just hardlinks to them from the docs
# directory.  This is just a little thing to make it easy to commit them into the
# docs directory for serving up on GitHub pages, while gitignoring them in the main directory.
#system("cd docs; rm [012]??*.html; for i in ../[012]??*.html; do echo $i; ln $i $(basename $i);  done")


