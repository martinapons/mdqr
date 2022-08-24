#' Create a result table for the regressor of interest
#'
#' @param object An object of class \code{mdqr}.
#' @param x A string containing the name of the regressor of interest. If left unspecified the regression results for all regressors are returned.
#' @return Table(s) of results.
#' @export

summary_mdqr <- function(object, x = NULL){
  if (is.null(x) == 1){
    coeff <- object[[1]]
  } else {
  res <- object[[1]]
  quantiles <- object[[3]]
  coeff <- matrix(NA, ncol = 4, nrow = length(quantiles))
  if (mode(res[[1]]) == "S4"){ # If GMM
    cnames <- c("Quantiles", colnames(res[[1]]@coef))
    for (u in quantiles){
    tab <- res[[which(u == quantiles)]]@coef # extract coefficient table
    coeff[which(u == quantiles),] <- tab[which(rownames(tab) == "x"),]
}
    } else if (methods::is(res[[1]],  "coeftest") == 1) { # for reoi
      cnames <- c("Quantiles", colnames(res[[1]]))
      for (u in quantiles){
        tab <- res[[which(u == quantiles)]] # extract coefficient table
        coeff[which(u == quantiles),] <- tab[which(rownames(tab) == "x"),]
      }
    } else {
   cnames <- c("Quantiles", colnames(res[[1]]$coeftable))
   for (u in quantiles){
     tab <- res[[which(u == quantiles)]]$coeftable # extract coefficient table
     coeff[which(u == quantiles),] <- tab[which(rownames(tab) == x),]
   }
  }

  coeff <- cbind(quantiles, coeff)
  colnames(coeff) <- cnames
  }
  return(coeff)
}
