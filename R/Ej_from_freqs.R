#' compute the expected posterior to a cluster from SNP allele freqs
#'
#' @param p vector of allele frequencies in K different populations
Ej_from_freqs <- function(p) {

  P00 <- p ^ 2
  P01 <- 2 * p * (1 - p)
  P11 <- (1 - p) ^ 2

  Q00 <- P00 / sum(P00)
  Q01 <- P01 / sum(P01)
  Q11 <- P11 / sum(P11)

  colSums(
    rbind(P00 * Q00,
          P01 * Q01,
          P11 * Q11
    )
  )
}
