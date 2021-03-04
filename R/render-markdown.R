# render an Rmarkdown document within an Rscript
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 1) {
  stop("Must be called from Rscript with path toe RMD file as the sole argument")
} else {
  f <- args[1]
}



outfile <- stringr::str_replace(f, "Rmd$", "stdout")
errfile <- stringr::str_replace(f, "Rmd$", "stderr")

message("")
message("Processing         : " , f)
message("Standard output to : ", outfile)
message("Error output to    : ", errfile)
message("")
# a bunch of rigamoral to get it to render in an entirely fresh, clean R session.
# See: https://github.com/rstudio/rmarkdown/issues/1204.  This is likely not
# necessary if running this within Snakemake.
path <- callr::r(
  function(...) rmarkdown::render(...),
  args = list(input = f,
              output_format = "html_document",
              envir = globalenv()
  ),
  stdout = outfile,
  stderr = errfile
)

