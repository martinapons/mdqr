# First stage function

md_first_stage <- function(data1 = subsample, time_varying_var, dep_var, U) {
  if (length(time_varying_var > 0)) { # if there are no time varying variables.
    ### form <- as.formula(paste(dep_var, paste(time_varying_var, collapse = "+"), sep = "~"))
    mm <- model.matrix(formula, data1) # make model matrix


    qr.X <- qr(mm, tol = 1e-9, LAPACK = FALSE) # remove collinear regressors
    rnkX <- qr.X$rank # (number of non-collinear columns)
    keepp <- qr.X$pivot[seq_len(rnkX)]
    mm <- mm[, keepp] # keep only non-collinear columns
    # form <- paste(dep_var, "~" , paste(colnames(mm)[-1], collapse  = "+"))
    # mm <- model.frame(formula(form), data1) # make model matrix
  } else {
    form <- as.formula(paste(dep_var, "1", sep = "~"))
    mm <- stats::model.matrix(formula(form), data1) # make model matrix
  }

  fitted <- array(NA, dim = c(dim(data1)[1], 1, length(U)))
  lambda_g <- array(0, dim = c(dim(mm)[2], dim(mm)[2], length(U))) # k = number of individual level regressors
  y <- model.frame(fdep, data1)[, 1]
  for (u in U) { # For all quantiles
    qreg <- rq(y ~ mm - 1, data = data1, tau = u)
    fitted[, 1, which(U == u)] <- fitted(qreg) # Save the fitted values
    lambda_g[, , which(U == u)] <- summary.rq(qreg, se = "ker", covariance = T)$cov / dim(data1)[1]
  }
  store <- list(fitted, lambda_g, mm)
  names(store) <- c("fitted", "lambda_g", "mm")
  return(store)
}
