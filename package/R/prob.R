#' @export
pvalOfBEZI <- function(x, mu, sigma, nu, lower = TRUE){
  pval <- (1 - gamlss.dist::pBEZI(
    x, mu = mu, sigma = sigma, nu = nu, lower.tail = lower, log.p = FALSE))
  return(pval)
}
