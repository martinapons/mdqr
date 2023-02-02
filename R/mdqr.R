#' Run minimum distance quantile regression.
#' @export
#'
#' @param formula An object of class \code{\link[Formula]{Formula}}. The formula consists of five parts.
#' \code{y ~ exo_1 + exo_2 | endo_1 + endo_2 | inst_1 + inst_2 | fe_1 + fe_2 | group_ID}.
#' Note that all part of the formula have to be specified. Use \code{0} to leave
#' it unspecified, e.g.\code{y ~ exo_1 + exo_2 | 0 | 0 | fe_1 + fe_2 | group_ID }.
#' Only second stage fixed effects should be included in the formula. If you wish to
#'  estimate a within model you don't need to include fixed effects (see Example 4).
#' @param data 	A data.frame containing the necessary variables to run the model.
#' @param method A character scalar indicates which method should be used in the second stage.
#'  The second stage estimators can use only the within variation (i.e. using an individual
#'   fixed effects approach) with the option "within", or the between variation with the option \code{"between"}.
#'    The options \code{"reoi"}, and \code{"regmm"} perform random effects estimation implemented with the optimal
#'     instrument and efficient GMM, respectively. \code{"ols"} performs a least-squares second stage of the
#'      fitted values on X. For \code{"gmm"}, the efficient weighting matrix is used. Second stage fixed effect
#'       might be included with \code{"ols"}, \code{"2sls"}, and \code{"gmm"}.
#' @param quantiles A vector with the quantiles of interest. The default is 0.1, 0.2, ... , 0.9.
#' @param clustervar A string with the name of the cluster variable. If \code{clustervar = NULL} (default),
#'  group indicator is used as a cluster variable.
#' @param cores Number of cores to use for first stage computation. if \code{core = NULL} the number of cores
#'  is set to 1. To see the number of cores available in your computer type \code{\link[parallel]{detectCores}}.
#' @param n_small A positive integer indicating the minimum size of groups allowed.
#'  Groups strictly smaller than \code{n_small - #first stage regression - 1} are dropped from the sample.
#'  If left unspecified it is set to 1.
#' @param run_second A logical evaluating to \code{TRUE} or \code{FALSE} indicating
#'  whether the second stage should be performed.
#' @param fitted_values A matrix containing the first stage fitted values.
#'  To use only if the function \code{mdqr} has been already run and only the second stage is different.
#'   For example, to change the clustering of the errors or to change the set of second-stage fixed effects.
#' @param run_time A logical evaluating to \code{TRUE} or \code{FALSE} indicating
#'  whether the computation time should be printed.
#' @details # Time-varying and time-constant variables / Individual-level and Group-level variables
#'
#' The formula automatically selects the regressors that have to be included in the first stage. If an endogenous variable is specified, either the within estimator is used or second-stage fixed effects have to be specified.
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#'
#' @examples
#'
#'
#' ## Example 1: Grouped data, no instrument no second stage fixed effects
#' # Generate artificial data.
#' rm(list = ls())
#' set.seed(1234)
#' # Generate 100 groups with 30 individuals each.
#' G <- 100
#' Ng <- 30
#' group <- rep(1:G, each = Ng)
#' # Generate a binary treatment variable.
#' treatment <- rep(rbinom(G, 1,0.5), each = Ng)
#' # Generate a group-level variable.
#' g_var <- rep(rnorm(G), each = Ng)
#' # Generate an individual-level variable.
#' i_var <- rnorm(Ng*G)
#' # Generate the outcome.
#' y <- 1 + treatment + g_var + i_var + rnorm(Ng * G) * (1 + 0.2 * treatment)
#' # Estimate the quantile treatment effect
#' fit <- mdqr(y~treatment+i_var+g_var | 0 | 0 | 0 | group, method = "ols",
#' data = data.frame(y, i_var, g_var, treatment, group), cores = 1)
#' # Compare with the true treatment effects
#' qnorm(seq(0.1,0.9,0.1)) * 0.2 + 1
#' # Now estimate the results at 19 quantiles and plot the treatment effects
#' fit <- mdqr(y~treatment+i_var+g_var | 0 | 0 | 0 | group, method = "ols",
#'  data = data.frame(y, i_var, g_var, treatment, group), quantiles = seq(0.05, 0.95, 0.05), cores = 1)
#' summary_mdqr(fit, "treatment")
#' plot_mdqr(fit, "treatment")
#'
#'
#'
#' # Example 2: Grouped data, instrumental variable, no second stage fixed effects
#' rm(list = ls())
#' set.seed(1234)
#' # Generate 100 groups with 30 individuals each.
#' G <- 200
#' Ng <- 30
#' group <- rep(1:G, each = Ng)
#'
#' # Generate the treatment, an instrumental variable, and (unobservable) group effects
#' sigma <- matrix(c(1, 0.5, 0, 0.5, 1, 0.5, 0, 0.5 ,1), ncol = 3)
#' mat <- MASS::mvrnorm(G, mu = c(0,0,0), Sigma = sigma)
#' g_effects <- rep(mat[,1], each = Ng)
#' treatment <- rep(mat[,2], each = Ng)
#' instrument <- rep(mat[,3], each = Ng)
#' # Generate a group-level variable.
#' g_var <- rep(rnorm(G), each = Ng)
#' # Generate an individual-level variable.
#' i_var <- rnorm(Ng*G)
#' # Generate the outcome.
#' y <- 1 + treatment + g_var + i_var + g_effects + rnorm(Ng * G) * (1 + 0.2 * treatment)
#' # Estimate the treatment effect without an instrumental variable: the result is inconsistent.
#' fit <- mdqr(y~treatment+i_var+g_var | 0 | 0 | 0 | group, method = "ols",
#' data = data.frame(y, i_var, g_var, treatment, group), cores = 1)
#' # Estimate the treatment effect with an instrumental variable.
#' fitIV <- mdqr(y~i_var+g_var | treatment | instrument | 0 | group, method = "2sls",
#'  data = data.frame(y, i_var, g_var, treatment,instrument, group), cores = 1)
#' summary_mdqr(fitIV, "fit_treatment")
#' # Compare with the true treatment effects
#' qnorm(seq(0.1,0.9,0.1)) * 0.2 + 1
#' # Plot results
#' plot_mdqr(fitIV, "fit_treatment")
#'
#'
#' # Example 3: Grouped data, no instrument, many fixed effects in the second stage
#' rm(list = ls())
#' set.seed(1234)
#'
#' # Generate data for 50 states, 20 years where groups are state x years cells.
#' # Each group has 30 observations.
#' S <- 50
#' Y <- 20
#' G <- S * Y
#' Ng <- 30
#' # Generate group identifier
#' group <- rep(1:G, each = Ng)
#' # Generate 50 states (groups) with 30 years of data.
#' state <- rep(1:S, each = Ng*Y)
#' state_fe <- rep(0.2 * rnorm(S), each = Ng*Y)
#' first_treated <- rep(round(runif(S) * 20), each = Ng*Y)
#' year <- rep(rep(1:Y, each = Ng), S)
#' year_fe <- year * 0.001
#' treated <- as.numeric(year >= first_treated)
#' state_char <- rep(rnorm(S*Y), each = Ng)
#' ind_char <- rnorm(Ng * G)
#' # Generate the outcome. Note that we have state and year fixed effects.
#' y <- 1 + treated + 0.5* state_char + 0.5 * ind_char + state_fe + year_fe +
#' rnorm(Ng*G)*(1 + treated * 0.2)
#' # Generate a data frame
#' dat <- data.frame(y, state, state_fe, year, year_fe, treated, state_char, ind_char, group)
#' # True quantile effect.
#' qnorm(seq(0.1,0.9,0.1)) * 0.2 + 1
#' # Estimate the treatment effects
#' fit <- mdqr(y~treated + state_char + ind_char  | 0 | 0 | state_fe + year_fe | group, method = "ols",
#'  data = dat, cores = 1)
#' # Cluster the standard errors at the state level
#' fit <- mdqr(y~treated + state_char + ind_char  | 0 | 0 | state_fe + year_fe | group, method = "ols",
#'  data = dat, clustervar = "state", cores = 1)
#' # Alternatively, without re-computing the first stage
#' fit <- mdqr(y~treated + state_char + ind_char  | 0 | 0 | state_fe + year_fe | group, method = "ols",
#'  data = dat, clustervar = "state",
#'  fitted_values = fit$data[, ncol(fit$data)- length(fit$quantiles):1+1], cores = 1)
#' # Result Table
#' summary_mdqr(fit, "treated")
#' # Plot Results
#' plot_mdqr(fit, "treated")
#'
#' # Example 4: Classical Panel Data
#' rm(list = ls())
#' set.seed(1234)
#' # Generate data for 100 individuals and 20 periods.
#' N <- 100
#' TT <- 20
#' id <- rep(1:100, each = TT)
#' # The individual effects are independent: re, fe and be are consistent.
#' alpha <- rep(rnorm(N), each = TT)
#' x <- rnorm(N*TT)
#' # Generate the outcome variable.
#' y <- x + alpha + rnorm(N*TT) * (1 + x * 0.2)
#' # True quantile treatment effects at the 0.1, 0.25, 0.5, 0.75 and 0.9 quantiles:
#' 1+qnorm(seq(0.1,0.9,0.1))*0.2
#' # Create dataset
#' dat <- data.frame(y, x, id)
#' # Compute the fixed-effects estimator.
#' fitfe <- mdqr(y ~ x | 0 | 0 | 0 | id, data = dat, method = "within", cores = 1)
#' # Results table.
#' summary_mdqr(fitfe, "fit_x")
#' # Plot results
#' plot_mdqr(fitfe, "fit_x")
#' # Compute the between estimator.
#' fitbe <- mdqr(y ~ x | 0 | 0 | 0 | id, data = dat, method = "be", cores = 1)
#' # Results table.
#' summary_mdqr(fitbe, "fit_x")
#' # Plot results.
#' plot_mdqr(fitbe, "fit_x")
#' # Compute the random effects estimator.
#' fitregmm <- mdqr(y ~ x | 0 | 0 | 0 | id, data = dat, method = "regmm", cores = 1)
#' # Results table.
#' summary_mdqr(fitregmm, "x")
#' # Plot results.
#' plot_mdqr(fitregmm, "x")
#' # Alternatively
#' fitreoi <- mdqr(y ~ x | 0 | 0 | 0 | id, data = dat, method = "reoi", cores = 1)
#' # Results table.
#' summary_mdqr(fitreoi, "x")
#' # Plot results.
#' plot_mdqr(fitreoi, "x")
#'
#'
#' @return
#' A list of four elements. The first element contains regression results for each quantile. The second element contains the matrix of fitted values from the first stage. The third element is the vector of quantiles. The fourth elements contains the number of groups.
#' see https://github.com/lrberge/fixest/blob/HEAD/R/ESTIMATION_FUNS.R for inspiration
#'
#' @author
#' Martina Pons, Blaise Melly
#'
#' @references \href{https://martinapons.github.io/files/MD.pdf}{Melly Blaise, Pons Martina (2022): "Minimum Distance Estimation of Quantile Panel Data Models"}.


mdqr <- function(formula,
                 data,
                 method = c( "ols", "within", "be", "reoi", "regmm", "2sls", "gmm"),
                 quantiles = seq(0.1, 0.9, 0.1),
                 clustervar = NULL,
                 cores = NULL,
                 n_small = NULL,
                 run_second = TRUE,
                 fitted_values = NULL,
                 run_time = FALSE
                 ) {

  method = match.arg(method) # check imput arguments

  if (is.null(cores) ==1 ) {
    message("The code is running on 1 core. If you want to use parallel computing specify the number of cores using the option 'cores'. To detect the number of cores on the current computer run 'detectCores()'")
  }

  if (is.data.frame(data) != TRUE) {
    stop("data must be a data.frame.")
  }

  if (is.null(n_small)) {
    n_small <- 1
  }

  start   <- Sys.time()
  formula <- Formula::as.Formula(formula)
  myvar   <- all.vars(formula)
  data <- dplyr::select(data, tidyselect::all_of(myvar), tidyselect::all_of(clustervar) ) # dataset containing only necessary columns.
  data %<>% dplyr::as_tibble() %>% tidyr::drop_na() # drop observations with missing values
  group <- stats::model.frame(formula(formula, lhs = 0, rhs = 5), data)
  group_id <- names(group) # groups are defined by this variable
  names(group) <- c("group")

  if  (group_id != "group") {
    if (sum(names(data) == "group") == 0 ) {
      data <- dplyr::bind_cols(data, group)
    } else {
      stop("One variable is named group. Change his name.")
    }
    data <- dplyr::select(data, -(group_id))
  }

  # recode the variable group such that is start at 1 and ends at G
  data %<>%
    dplyr::as_tibble() %>%
    dplyr::group_by_at(dplyr::vars(group)) %>%  # Need to change this line! End of life cycle
    dplyr::mutate(group = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  data <- data %>% dplyr::add_count(group) # number of observations in each group.
  data <- data %>% dplyr::filter(n >= n_small + 1) # remove groups with less than n_small obs + 1 (at least 1 degree of freedom for the constant)

  # recode the variable group such that is start at 1 and ends at G
  data %<>%
    dplyr::as_tibble() %>%
    dplyr::group_by_at(dplyr::vars(group)) %>%
    dplyr::mutate(group = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  data <- dplyr::arrange(data, group) # sort the data
  G <- max(data$group) # number of groups

  group <- data$group
  if (is.null(clustervar)){ # if no clustervar is specified, the group identifies is used.
    clvar <- group
  } else {
    clvar <- dplyr::select(data, tidyselect::all_of(clustervar))
  }
  # ------------------------------------------------------------
  # Chunks of formula
  fdep <- formula(formula, lhs = 1, rhs = 0)
  fex  <- formula(formula, lhs = 0, rhs = 1)
  fen  <- formula(formula, lhs = 0, rhs = 2)
  ffe  <- formula(formula, lhs = 0, rhs = 4)
  fz   <- formula(formula, lhs = 0, rhs = 3)

  # Generate matrices with dependent variables
  y   <- stats::model.frame(fdep, data) # some of these are not needed
  end <- stats::model.frame(fen, data)
  z   <- stats::model.frame(fz, data)
  fe  <- stats::model.frame(ffe, data)
  exo <- stats::model.frame(fex, data)

  # -----------------

  endog_var <- names(end)
  exo_var   <- names(exo)

  if (length(all.vars(fen)) > length(all.vars(fz)) & method != "within" & method != "ht") stop("Fewer instruments than endogenous variables. If you wish to use interval instrument select method fe or ht")
  if (length(all.vars(fen)) == 0 & length(all.vars(fen)) > 0) stop("External instrument is specified, but there is no endogeous variable")
  if (method == "ols" & (length(all.vars(fen) ) >  0  & length(all.vars(ffe)) == 0 )) stop("OLS is used but there are endogenous variables.")
  if (method == "reoi" & length(all.vars(ffe)) > 0) stop("RE cannot be uesd with fixed effects in the second stage.")
  if (method == "ht" & length(all.vars(ffe)) > 0) stop("HT cannot be used with fixed effects in the second stage.")
  if (method == "regmm" & length(all.vars(ffe)) > 0) stop("RE cannot be uesd with fixed effects in the second stage.")
  if (method == "within" & length(all.vars(ffe)) > 0) stop("The within estimator cannot be use with fixed effects in the second stage. The within estimator only exploit variation within individuals (groups). If you want to include fixed effects as a higher level than the individual (group), use the option ols, gmm, or iv")
  if (length(all.vars(fz)) > 0){
    for (l in 1:length(all.vars(fz))){
      if (sum(tapply(z[[l]],data$group, stats::var) != 0) != 0 )  stop("The instrument is varying within individuals (groups). The instrument is only allowed to vary between individuals (groups)")
    }
  }
  if (run_second == F & (is.null(fitted_values) == 0)) stop("No estimation to perform. Either you set run_second = TRUE or you let the fitted values unspecified.")
  # --------------------------------------------------------------
  if (is.null(fitted_values) == 1) {
    if (run_time == TRUE){
    }
    # This part of the code generates a list containing the data for each first stage regression.
    datalist <- data %>% dplyr::group_split(group)

    # First stage ---------------------------------------------------------------------------------
    if (run_time == TRUE){
      message("First stage estimation starting...")
      print(Sys.time()- start)

    }
    if (is.null(cores) == 1) {
      cores <- 1
    }
    form1 <- formula(paste0(as.character(fdep)[2], "~", as.character(fex)[2], "+ ", (as.character(fen)[2])))


      future::plan(future::multisession, workers = cores, gc = T)

    #options(future.globals.maxSize = 8000 * 1024^2)
    first <- furrr::future_map(datalist, md_first_stage, fdep, form1, quantiles, n_small
    )
    # cl <- parallel::makeCluster(cores) # Set the number of clusters
    # parallel::clusterExport(cl, c("fdep", "md_first_stage", "form1", "datalist", "quantiles", "n_small"), envir=environment()) # Functions needed
    # parallel::clusterEvalQ(cl, {
    #   library(quantreg)
    # })
    # environment(fdep) <- .GlobalEnv
    # environment(form1) <- .GlobalEnv
    # environment(quantiles) <- .GlobalEnv
    # environment(md_first_stage) <- .GlobalEnv
    # environment(datalist) <- .GlobalEnv
    # environment(n_small) <- .GlobalEnv
    #
    # first <- parallel::clusterApply(cl, datalist, md_first_stage, fdep, form1, quantiles, n_small)
    #
    # parallel::stopCluster(cl)

    # Extract Fitted values: the fitted values are a list in a list.
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

  if (run_second == TRUE) {
    if (run_time == TRUE){
      message("Second stage estimation starting...")
      print(Sys.time()- start)
    }
    if (fex == "~0") {
      fex <- "~1"
    }
    if (method == "within") { # no second stage FE possible method = "within"
      b <- cbind(group, exo) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~ . - mean(.)) , .names = "{.col}dem" ))
      inst_s <- paste0(paste(exo_var, collapse = "dem +"), "dem")
      second <- dplyr::bind_cols(b[, -1], data, fitted)
      form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~ 1" ,  "|", as.character(fex)[2], "~", inst_s))

      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "ols") {
      if (length(fe) == 0) {
        form <- Formula::as.Formula(paste0(".[mydep]~", as.character(fex)[2]))
      } else {
        form <- Formula::as.Formula(paste0(".[mydep]", paste0("~", as.character(fex)[2], "|", as.character(ffe)[2])))

      }
      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "2sls") {
      if (length(fe) == 0) {
        form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~" , as.character(fex)[2], "|", as.character(fen)[2], "~", as.character(fz)[2]))
      } else {
        form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~" , as.character(fex)[2], "|", as.character(ffe)[2], "|", as.character(fen)[2], "+", as.character(fz)[2]))

      }
      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "be") {
      b <- cbind(group, exo) %>%
        dplyr::group_by(group) %>%
        dplyr::transmute(dplyr::across(tidyselect::everything(),  list(tdm = ~ mean(.)) , .names = "{.col}m" ))
      second <- dplyr::bind_cols(b[, -1], data, fitted)
      inst_s <- paste0(paste(exo_var, collapse = "m +"), "m")
      form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~1", "|", as.character(fex)[2], "~", inst_s))

      res <-  fixest::feols(form, second, cluster = clvar)
    } else if (method == "reoi") {
      lambda <- sapply(first, function(x) x[c("lambda_i")])
      xlist <- sapply(first, function(x) x[c("mm")])
      Xtilde <- Matrix::bdiag(xlist)
      res <- list()
      for (u in quantiles) {
        Lambda <- lapply(lambda, function(x) x[, , which(u == quantiles)])
        # this chunk of code comes from https://github.com/cran/plm/blob/master/R/tool_ercomp.R. This uses Nerlove method and
        # with OLS using the fitted values lead to the same estimated sigma_alpha.

        form <- Formula::as.Formula(paste0(paste0("fitted_", u , "~"), as.character(fex)[2]))
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
        form2 <- Formula::as.Formula(paste0(paste0("fitted_", u , "~"), as.character(fex)[2], "| Zstar-1 "))
        rr <- AER::ivreg(form2, data = second)
        rr <- lmtest::coeftest(rr, vcov = sandwich::vcovCL, cluster = clvar)
        res[[which(u == quantiles)]] <- rr
      }
      names(res) <- mydep
    } else if (method == "gmm") {
      # https://cran.r-project.org/web/packages/momentfit/vignettes/gmmS4.pdf
      res <- list()

      for (u in quantiles) {
        form <- Formula::as.Formula(paste0(paste0("fitted_", u), "~" , as.character(fex)[2] , "+", as.character(fen)[2]))

        model <- momentfit::momentModel(form, fz, data = second, vcov = "CL", vcovOptions = list(cluster = data.frame(clvar)))
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

      Z <- cbind(1, de, me)
      res <- list()
      for (u in quantiles) {
        form <- Formula::as.Formula(paste0(paste0("fitted_", u),  "~" , as.character(fex)[2]))
        model <- momentfit::momentModel(form, Z, data = second, vcov = "CL", vcovOptions = list(cluster = data.frame(clvar)))
        r <- momentfit::gmmFit(model, type = "twostep")
        rr <- momentfit::summary(r, sandwich = TRUE, df.adj = FALSE) # see 1.4.4 in https://cran.r-project.org/web/packages/momentfit/vignettes/gmmS4.pdf
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
      form <- Formula::as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", as.character(fex)[2], "|",  as.character(fen)[2], inst_s))

      res <-  fixest::feols(form, second, cluster = clvar)
    }
  } else {
    res <- NULL
  }
  # count group that were used for estimation
  G <- length(unique(group[!is.na(fitted[, 1])]))

  # save results into a list.
  res <- list(res, second, quantiles, G)
  names(res) <- c("results", "data", "quantiles", "G" )
  if (run_time == TRUE){
    print(Sys.time()- start)
  }

  class(res) <- "mdqr"


  return(res)

}


