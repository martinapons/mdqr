# First stage function

md_first_stage <- function(data1 = subsample, fdep, form1, U) {

    mm <- model.matrix(form1, data1) # make model matrix
    qr.X <- qr(mm, tol = 1e-9, LAPACK = FALSE) # remove collinear regressors
    rnkX <- qr.X$rank # (number of non-collinear columns)
    keepp <- qr.X$pivot[seq_len(rnkX)]
    mm <- mm[, keepp] # keep only non-collinear columns


  fitted <- array(NA, dim = c(dim(data1)[1], 1, length(U)))

  lambda_i <- array(0, dim = c(dim(mm)[2], dim(mm)[2], length(U))) # k = number of individual level regressors
  if (is.null(c(dim(mm)[2]))){
    lambda_i <- array(0, c(1,1, length(U)))
  }
  y <- model.frame(fdep, data1)[, 1]
  for (u in U) { # For all quantiles
    qreg <- rq(y ~ mm - 1, data = data1, tau = u)
    fitted[, 1, which(U == u)] <- fitted(qreg) # Save the fitted values
    lambda_i[, , which(U == u)] <- summary.rq(qreg, se = "ker", covariance = T)$cov / dim(data1)[1]
  }
  store <- list(fitted, lambda_i, mm)
  names(store) <- c("fitted", "lambda_i", "mm")
  return(store)
}
