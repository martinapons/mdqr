#' Run first stage regression
#' This function is called by the mdqr() function.
#' @export
#' @param data1 data Subsample of dataset containing one group.
#' @param fdep formula containing the dependent variable.
#' @param form1 formula containing endogenous and exogenous regressors.
#' @param quantiles vector of quantiles.

md_first_stage <- function(data1, fdep, form1, quantiles) {

  m <- stats::model.frame(form1, data1) # this deals with cases with factor variables with 1 level only.
  # esclude variables that do not vary within groups.
  var <- apply(m, 2, var)
  sel <- names(var[var != 0])
  if (var[1] == 0){ # Check that there is within-group variation in the dependent variable.
    stop("there is no variation in the dependent variable in at least one group")
  }
  if (length(sel) == 1){ # if only the dependent variable varies within group, the first stage consist of a quantile regression on a constant.
    form1 <- stats::as.formula(paste0(sel[1], "~",1 ))
    mm <- stats::model.matrix(form1, data1)
  } else {
    form1 <- stats::as.formula(paste0(sel[1], "~", paste(sel[2:(length(sel))], collapse = "+")))

    mm <- stats::model.matrix(form1, data1) # create model matrix

    # This code check that there are no collinear regressors. This could happen, e.g. with some factor variables or some interaction.
    qr.X <- qr(mm, tol = 1e-9, LAPACK = FALSE)
    rnkX <- qr.X$rank
    keepp <- qr.X$pivot[seq_len(rnkX)]
    mm <- mm[, keepp]
  }

  fitted <- array(NA, dim = c(dim(data1)[1], 1, length(quantiles))) # Create array to store the fitted values

  lambda_i <- array(0, dim = c(dim(mm)[2], dim(mm)[2], length(quantiles))) # Create array to store the first stage variances. k = number of individual level regressors
  if (is.null(c(dim(mm)[2]))){
    lambda_i <- array(0, c(1,1, length(quantiles)))
  }
  y <- stats::model.frame(fdep, data1)[, 1]
  for (u in quantiles) {
    qreg <- quantreg::rq(y ~ mm - 1, data = data1, tau = u) # first stage quantile regression
    fitted[, 1, which(quantiles == u)] <- fitted(qreg) # Save the fitted values
    lambda_i[, , which(quantiles == u)] <- quantreg::summary.rq(qreg, se = "ker", covariance = T)$cov / dim(data1)[1] # covariance matrix for each group
  }
  store <- list(fitted, lambda_i, mm)
  names(store) <- c("fitted", "lambda_i", "mm")
  return(store)
}
