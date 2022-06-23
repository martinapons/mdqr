#' Create a result table for the regressor of interest
#' @export
#'
#' @param object An object estimated by \code{mdqr}.
#' @param x A string containing the name of the regressor of interest.
#' @return Table of results.

result_table <- function(object, x){
  res <- object[[1]]
  quantiles <- object[[3]]
  coeff <- matrix(NA, ncol = 4, nrow = length(quantiles))
  cnames <- c("Quantiles", colnames(res[[1]]$coeftable))
  for (u in quantiles){
    tab <- res[[which(u == quantiles)]]$coeftable
    coeff[which(u == quantiles),] <- tab[which(rownames(tab) == x),]
  }
  coeff <- cbind(quantiles, coeff)
  colnames(coeff) <- cnames
  return(coeff)
}
