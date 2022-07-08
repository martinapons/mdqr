library(parallel)
library(stringr)
library(Formula)
library(MASS)
library(emulator)
library(dplyr)
library(fixest)
library(tidyselect)
library(tidyr)
library(Matrix)
quantile <- c(.1, .5)

data <- d

#' Run minimum distance regression.
#' @export
#'
#' @param formula An object of class \code{\link[Formula]{Formula}}. The formula consists of five parts. \code{y ~ exo_1 + exo_2 | endo_1 + endo_2 | z_1 + z_2 | fe_1 + fe_2 | group_ID }. Note that all part of the formula have to be specified. Use \code{0} to leave it unspecified, e.g.\code{y ~ exo_1 + exo_2 | 0 | 0 | fe_1 + fe_2 | group_ID }. Only second stage fixed effects should be included in the formula. If you wish to estimate a within model you don't need to include fixed effects (see 'details').
#' @param data 	A data.frame containing the necessary variables to run the model.
#' @param method A character scalar indicates which method should be used in the second stage. The second stage estimators can use only the within variation (i.e. using an individual fixed effects approach) "within", the between variation "between". "reoi", and "regmm" perform random effects estimation implemented with the optimal instrument and efficient GMM, respectively. "ols" performs a least-squares second stage of the fitted values on X. "reoi" and "regmm" perform a second stage regression of the fitted values on X using an instrument vector Z. For "gmm", the efficient weighting matrix is used. Second stage fixed effect might be included with "ols", "2sls", and "gmm".
#' @param quantiles A vector with the quantiles of interest. The default is 0.1, 0.2, ... , 0.9.
#' @param clustervar A string with the name of the cluster variable. If \code{clustervar = NULL} (default), group indicator is used as a cluster variable.
#' @param cores Number of cores to use for first stage computation. if \code{core = NULL} the number of cores is set to \code{\link[parallel]{detectCores}-1}.
#' @param n_small A positive integer indicating the minimum size of groups allowed. Groups strictly smaller than \code{n_small} are dropped from the sample.
#' @param run_second A logical evaluating to \code{TRUE} or \code{FALSE} indicating whether the second stage should be performed.
#' @param fitted_values A matrix containing the first stage fitted values. To use only if the function \code{mdqr} has been already run and only the second stage is different. For example, to change the clustering of the errors or to change the set of second stage fixed effects.
#' @section Time-varying and time-constant variables / Individual-level and Group-level variables
#' The formula automatically select the regressors that have to be included in the first stage. If an endogenous variables is specified, either the within estimator is used or second stage fixed effects have to be specified.
#' @section Implementing the different estimators

#' ## Within Regression:
#'
#' Estimate the effect of union status on wage using a fixed effects regression.
#' \code{formula(wage ~ union | 0 | 0 | 0 | group_ID, data, method = c("within"))}
#'
#' or
#'
#' #' \code{formula(wage ~ 0 | union | 0 | 0 | group_ID, data, method = c("within"))}
#'
#' ## Random Effects Regression
#'
#' \code{formula(wage ~ 0 | union | 0 | 0 | group_ID, data, method = c("reoi"))}
#' \code{formula(wage ~ 0 | union | 0 | 0 | group_ID, data, method = c("regmm"))}
#'
#' ## Instrumental Variables
#' Assume we want to estimate the effect of school budget (x) students' outcomes (y), where schools define the groups. Assume that an instrument z is available.
#' \code{formula(y ~ 0 | x | z | 0 | group_ID, data, method = c("2sls"))}
#' Covariates and or say county fixed effects can be included as follows:
#' \code{formula(y ~ w | x | z | county | group_ID, data, method = c("2sls"))}
#'
#' If the model is overidentified "gmm" can be used instead.
#'
#' @return
#' A list of four elements. The first element contains regression results for each quantile. The second element contains the matrix of fitted values from the first stage. The third element is the vector of quantiles.
#' see https://github.com/lrberge/fixest/blob/HEAD/R/ESTIMATION_FUNS.R for inspiration
#' @author
#' Martina Pons
#'
#' @references Melly Blaise, Pons Martina (2022): "Minimum Distance Estimation of Quantile Panel Data Models" ([](https://martinapons.github.io/files/MD.pdf)).
mdqr <- function(formula, data, method = c("within", "be", "reoi", "regmm", "ols", "2sls", "gmm"), quantiles = seq(0.1, 0.9, 0.1), clustervar = NULL, cores = NULL, n_small = 1, run_second = TRUE, fitted_values = NULL) {
  start <- Sys.time()
  formula <- Formula::as.Formula(formula)

  myvar <- all.vars(formula)
  data <- dplyr::select(data, tidyselect::all_of(myvar), tidyselect::all_of(clustervar) ) # eventually need to add cluster var.

  data %<>% dplyr::as_tibble() %>% tidyr::drop_na()

  group <- stats::model.frame(formula(formula, lhs = 0, rhs = 5), data)

  group_id <- names(group) # groups are define by this variable

  names(group) <- c("group")


  # drop too small groups ----------------------------------------

  if  (group_id != "group") {
    if (sum(names(data) == "group") == 0 ) {
      data <- dplyr::bind_cols(data, group)
    } else {
      stop("One variable is named group. Change his name.")
    }
    data <- dplyr::select(data, -(group_id))
  }

  data %<>%
    dplyr::as_tibble() %>%
    dplyr::group_by_at(dplyr::vars(group)) %>%  # Need to change this line! End of life cycle
    dplyr::mutate(group = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  data <- data %>% dplyr::add_count(group) # number of observations in each group.
  data <- data %>% dplyr::filter(n >= n_small) # remove groups with less than n obs

  data %<>%
    dplyr::as_tibble() %>%
    dplyr::group_by_at(dplyr::vars(group)) %>%
    dplyr::mutate(group = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  data <- dplyr::arrange(data, group)
  G <- max(data$group)

  if (is.null(clustervar)){
    clvar <- group
  } else {
    clvar <- dplyr::select(data, tidyselect::all_of(clustervar))
  }
  # ------------------------------------------------------------

  fdep <- formula(formula, lhs = 1, rhs = 0)
  fex <- formula(formula, lhs = 0, rhs = 1)
  fen <- formula(formula, lhs = 0, rhs = 2)
  ffe <- formula(formula, lhs = 0, rhs = 4)
  fz <- formula(formula, lhs = 0, rhs = 3)

  y <- stats::model.frame(fdep, data) # some of these are not needed
  end <- stats::model.frame(fen, data)
  z <- stats::model.frame(fz, data)
  fe <- stats::model.frame(ffe, data)
  exo <- stats::model.frame(fex, data)

  # -----------------

  endog_var <- names(end)
  exo_var <- names(exo)


  if (length(all.vars(fen)) > length(all.vars(fen)) & method != "within" & method != "ht") stop("Fewer instruments than endogenous variables. If you wish to use interval instrument select method fe or ht")
  if (length(all.vars(fen)) == 0 & length(all.vars(fen)) > 0) stop("External instrument is specified, but there is no endogeous variable")
  # if (min(tapply(y[, 1], group[, 1], var) > 0) == 0) stop("The dependent variable must vary within groups / individuals.")
#  if (length(all.vars(fen)) == 0 & method == "within") stop("The within estimator is used, but no endogenous variable specified.")
  if (method == "ols" & (length(all.vars(fen) ) >  0  & length(all.vars(ffe)) == 0 )) stop("OLS is used but there are endogenous variables.")
  if (method == "reoi" & length(all.vars(ffe)) > 0) stop("RE cannot be uesd with fixed effects in the second stage.")
  if (method == "ht" & length(all.vars(ffe)) > 0) stop("HT cannot be used with fixed effects in the second stage.")
  if (method == "regmm" & length(all.vars(ffe)) > 0) stop("RE cannot be uesd with fixed effects in the second stage.")
  if (method == "within" & length(all.vars(ffe)) > 0) stop("The within estimator cannot be use with fixed effects in the second stage. The within estimator only exploit variation within individuals (groups). If you want to include fixed effects as a higher level than the individual (group), use the option ols, gmm, or iv")

  if (length(all.vars(fz)) > 0){
    if (sum(tapply(z[[1]], group, stats::var) != 0) != 0 )  stop("The instrument is varying within individuals (groups). The instrument is only allowed to vary between individuals (groups)")
  }

  # --------------------------------------------------------------
  if (is.null(fitted_values) == 1) {
  print(Sys.time()-start)

  print("Make datalist...")

  first <- data
  datalist <- NULL
  first$g <- round(first$group * 0.001) + 1
  gg <- max(first$g)

  for (i in 1:gg) {
    first1 <- first %>% dplyr::filter(g == i)
    gr <- c(min(first1$group):max(first1$group))
    datalist1 <- lapply(gr, function(gr) first1 %>% dplyr::filter(group == gr))
    datalist <- c(datalist, datalist1)
  }

  # First stage ---------------------------------------------------------------------------------
  print(Sys.time()- start)

  print("First stage estimation starting...")
  if (is.null(cores) == 1) {
    cores <- parallel::detectCores() - 1
  }

  #form1 <- formula(paste0(stringr::str_sub(tchar(fdep), 1, -5), paste0(fex), "+ ", stringr::str_sub(tchar(fen), 2, -1)))
  form1 <- formula(paste0(as.character(fdep)[2], "~", as.character(fex)[2], "+ ", (as.character(fen)[2])))

  cl <- parallel::makeCluster(cores) # Set the number of clusters
  parallel::clusterExport(cl, c("fdep", "md_first_stage", "form1", "quantiles"), envir=environment()) # Functions needed
  parallel::clusterEvalQ(cl, {
    library(quantreg)
  })
  environment(fdep) <- .GlobalEnv
  environment(form1) <- .GlobalEnv
  environment(quantiles) <- .GlobalEnv
  environment(md_first_stage) <- .GlobalEnv
  environment(datalist) <- .GlobalEnv


  print("really starting...")
  first <- parallel::clusterApply(cl, datalist, md_first_stage, fdep, form1, quantiles)

  parallel::stopCluster(cl)
  print("parallel done...")

  print(Sys.time()- start)

  # Extract Fitted values: the fitted values are a list in a list. E
  fitted <- lapply(first, function(x) x[c("fitted")])
  fitted <- lapply(fitted, data.frame, stringsAsFactors = FALSE)
  fitted <- dplyr::bind_rows(fitted)
  # -------------------

  mydep <- paste("fitted", quantiles, sep = "_")
  colnames(fitted) <- mydep
  } else {
    fitted <- fitted_values
    mydep <- paste("fitted", quantiles, sep = "_")
    if (sum(colnames(fitted) != mydep) != 0) stop("The matrix of fitted values was either not generated by mdqr or was estimated for different quantiles.")
  }
  second <- dplyr::bind_cols(data, fitted)
  print(Sys.time()- start)

  if (run_second == TRUE) {
    print("Second stage estimation starting...")


    if (fex == "~0") {
      fex <- "~1"
    }
    if (method == "within") { # no second stage FE possible if you have  method = "within"
      b <- cbind(group, end) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~ . - mean(.)) , .names = "{.col}dem" ))
      inst_s <- paste0(paste(endog_var, collapse = "dem +"), "dem")
      second <- dplyr::bind_cols(b[, -1], data, fitted)
      #form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", stringr::str_sub(tchar(fen), 2, -1), "~", inst_s))
      form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~" , as.character(ffe)[2], "|", as.character(fen)[2], "~", inst_s))
      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "ols") {
      if (length(fe) == 0) {
        form <- Formula::as.Formula(paste0(".[mydep]", fex))
      } else {
        #form <- Formula::as.Formula(paste0(".[mydep]", paste(fex, "|", stringr::str_sub(tchar(ffe), 2, -1))))
        form <- Formula::as.Formula(paste0(".[mydep]", paste0("~", as.character(fex)[2], "|", as.character(ffe)[2])))
      }
      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "2sls") {
      if (length(fe) == 0) {
        #form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", stringr::str_sub(tchar(fen), 2, -1), fz))
        form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~" , as.character(fex)[2], "|", as.character(fen)[2], fz))
      } else {
        #form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", stringr::str_sub(tchar(ffe), 2, -1), "|", stringr::str_sub(tchar(fen), 2, -1), fz))
        form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~" , as.character(fex)[2], "|", as.character(ffe)[2], "|", as.character(fen)[2], fz))
      }
      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "be") {
      b <- cbind(group, end) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~ mean(.)) , .names = "{.col}m" ))
      second <- dplyr::bind_cols(b[, -1], data, fitted)
      inst_s <- paste0(paste(endog_var, collapse = "m +"), "m")
      #form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", stringr::str_sub(tchar(fen), 2, -1), "~", inst_s))
      form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~" , as.character(fex)[2], "|", as.character(fen)[2], "~", inst_s))
      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "reoi") {
      lambda <- sapply(first, function(x) x[c("lambda_i")])
      xlist <- sapply(first, function(x) x[c("mm")])
      Xtilde <- Matrix::bdiag(xlist)
      res <- list()
      for (u in quantiles) {
        Lambda <- lapply(lambda, function(x) x[, , which(u == quantiles)])
        # this chunk of code comes from https://github.com/cran/plm/blob/master/R/tool_ercomp.R. This uses Nerlove method and
        # with OLS using the fitted values lead to the same estiamted sigma_alpha.

        form <- Formula::as.Formula(paste0(paste0("fitted_", u), fex))
        est <- plm::plm(form, data = data.frame(second), index = c("group"), model = "within")
        pdim <- plm::pdim(est)
        N <- pdim$nT$n

        s2alpha <- sum((plm::fixef(est, type = "dmean", effect = "individual"))^2 *
                         pdim$Tint$Ti / pdim$nT$N) * (pdim$nT$n / (pdim$nT$n - 1))

        xlx <- Map(emulator::quad.tform, Lambda, xlist)
        omega <- Map("+", xlx, s2alpha)
        omegainv <- lapply(omega, MASS::ginv)
        Omegainv <- Matrix::bdiag(omegainv)
        Zstar <- as.matrix(Omegainv %*% stats::model.matrix(fex, data))
        form2 <- Formula::as.Formula(paste0(form, "| Zstar-1 "))
        rr <- AER::ivreg(form2, data = second)
        rr <- lmtest::coeftest(rr, vcov = sandwich::vcovCL, cluster = clvar)
        res[[which(u == quantiles)]] <- rr
      }
      names(res) <- mydep
    } else if (method == "gmm") {
      # https://cran.r-project.org/web/packages/momentfit/vignettes/gmmS4.pdf
      res <- list()

      for (u in quantiles) {
        #form <- Formula::as.Formula(paste0(paste0("fitted_", u), fex, "+", stringr::str_sub(tchar(fen), 2, -1)))

        form <- Formula::as.Formula(paste0(paste0("fitted_", u), "~" , as.character(fex)[2] , "+", as.character(fen)[2]))
        model <- momentfit::momentModel(form, fz, data = second, vcov = "CL", vcovOptions = list(cluster = clvar))
        rr <- summary(momentfit::gmmFit(model, type = "twostep"))
        res[[which(u == quantiles)]] <- rr
      }
      names(res) <- mydep
    } else if (method == "regmm") {
      me <- cbind(group, exo) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~ mean(.)) , .names = "{.col}m" ))

      dm <- cbind(group, exo) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~ . - mean(.)) , .names = "{.col}dem" ))

      me <- me[, -1]
      de <- dm[, -1]

      # drop zero in demeaned
      tv <- (colSums((abs(de) == 0)) == 0)
      de <- de[tv]

      Z <- cbind(de, me)

      for (u in quantiles) {
        #form <- Formula::as.Formula(paste0(paste0("fitted_", u), fex, "+", stringr::str_sub(tchar(fen), 2, -1)))
        form <- Formula::as.Formula(paste0(paste0("fitted_", u),  "~" , as.character(fex)[2], "+",  as.character(fen)[2]))

        model <- momentfit::momentModel(form, Z, data = second, vcov = "CL", vcovOptions = list(cluster = clvar))
        r <- momentfit::gmmFit(model, type = "twostep")
        rr <- summary(r, sandwich = TRUE, df.adj = FALSE) # whichone do we want? see 1.4.4 in https://cran.r-project.org/web/packages/momentfit/vignettes/gmmS4.pdf
        res[[which(u == quantiles)]] <- rr
      }
      names(res) <- mydep
    } else if (method == "ht") {
      me <- cbind(group, exo) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~  mean(.)) , .names = "{.col}m" ))

      dm <- cbind(group, exo, end) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~ . - mean(.)) , .names = "{.col}dem" ))
      me <- me[, -1]
      de <- dm[, -1]

      # drop zero in demeaned
      tv <- (colSums((abs(de) == 0)) == 0)
      de <- de[tv]
      Z <- cbind(de, me)
      second <- cbind(second, Z)
      inst_s <- paste0( "~" ,paste(names(Z), collapse = "+"))
      #form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", stringr::str_sub(tchar(fen), 2, -1), inst_s))
      form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", s.character(fex)[2], "|",  as.character(fen)[2], inst_s))
      res <-  fixest::feols(form, second, cluster = clvar)
    }
  } else {
    res <- NULL
  }
  res <- list(res, fitted, quantiles, G)
  names(res) <- c("results", "fitted_values", "quantiles", "G" )
  print("total time: ")
  print(Sys.time()- start)
  return(res)
  print(res[[1]])


}


