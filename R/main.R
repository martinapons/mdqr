library(purrrlyr)
library(parallel)
library(stringr)
library(Formula)
library(formula.tools)
library(swCRTdesign)
library(MASS)
library(emulator)
library(dplyr)
library(fixest)
library(tidyverse)
library(Matrix)
quantile <- c(.1, .5)

data <- d
mdqr <- function(formula, data, method = c("fe", "be", "reoi", "regmm", "ols", "2sls", "gmm"), quantiles = seq(0.1, 0.9, 0.1), cores = NULL) {
  formula <- Formula::as.Formula(formula)

  myvar <- all.vars(formula)
  data <- dplyr::select(data, tidyselect::all_of(myvar)) # eventually need to add cluster var.
  data %<>% as_tibble() %>% drop_na()

  group <- model.frame(formula(formula, lhs = 0, rhs = 5), data)
  group_id <- names(group) # groups are define by this variable
  names(group) <- c("group")

  # drop too small groups ----------------------------------------
  data %<>%
    as_tibble() %>%
    group_by_at(vars(group)) %>%
    mutate(group = cur_group_id()) %>%
    ungroup()

  data <- data %>% add_count(group) # number of observations in each group.
  data <- data %>% filter(n >= 2) # remove groups with less than 25 obs

  data %<>%
    as_tibble() %>%
    group_by_at(vars(group)) %>%
    mutate(group = cur_group_id()) %>%
    ungroup()

  G <- max(data$group)
  # ------------------------------------------------------------

  fdep <- formula(formula, lhs = 1, rhs = 0)
  fex <- formula(formula, lhs = 0, rhs = 1)
  fen <- formula(formula, lhs = 0, rhs = 2)
  ffe <- formula(formula, lhs = 0, rhs = 4)
  fz <- formula(formula, lhs = 0, rhs = 3)

  y <- model.frame(fdep, data) # some of these are not needed
  end <- model.frame(fen, data)
  z <- model.frame(fz, data)
  fe <- model.frame(ffe, data)
  exo <- model.frame(fex, data)

  # -----------------

  endog_var <- names(end)
  exo_var <- names(exo)


  if (length(all.vars(fen)) > length(all.vars(fen)) & method != "fe" & method != "ht") stop("fewer instruments than endogenous variables. If you wish to use interval instrument select method fe or ht")
  if (length(all.vars(fen)) == 0 & length(all.vars(fen)) > 0) stop("External instrument is specified, but there is no endogeous variable")
  if (min(tapply(y[, 1], group[, 1], var) > 0) == 0) stop("The dependent variable must vary within groups / individuals.")
  if (length(all.vars(fen)) == 0 & method == "fe") stop("FE is used, but no endogenous variable specified.")
  if (method == "ols" & length(all.vars(fen)) > 0) stop("OLS is used but there are endogenous variables.")
  if (method == "reoi" & length(all.vars(ffe)) > 0) stop("RE is uesd with fixed effects in the second stage.")

  if (sum(cbind(z, group) %>% slice_rows("group") %>% dmap(var) %>% colSums() != 0) != 1) {
    stop("The instrument is not allowed to vary within groups.")
  }
  # --------------------------------------------------------------


  first <- data
  datalist <- NULL
  first$g <- round(first$group * 0.001) + 1
  gg <- max(first$g)

  for (i in 1:gg) {
    first1 <- first %>% filter(g == i)
    gr <- c(min(first1$group):max(first1$group))
    datalist1 <- lapply(gr, function(gr) first1 %>% filter(group == gr))
    datalist <- c(datalist, datalist1)
  }

  # First stage ---------------------------------------------------------------------------------
  if (is.null(cores) == 1) {
    cores <- detectCores() - 1
  }

  form1 <- formula(paste0(stringr::str_sub(as.character(fdep), 1, -5), paste0(fex), "+ ", stringr::str_sub(as.character(fen), 2, -1)))

  cl <- makeCluster(cores) # Set the number of clusters
  RNGkind(kind = "L'Ecuyer-CMRG") # Set the seed for each cluster
  clusterExport(cl, c("fdep", "md_first_stage", "datalist", "form1", "U")) # Functions needed
  clusterEvalQ(cl, {
    library(quantreg)
    library(sandwich)
    library(dtplyr)
    library(pracma)
    library(stringr)
    library(Formula)
    library(Matrix)
  })
  first <- clusterApply(cl, datalist, md_first_stage, fdep, form1, U)
  stopCluster(cl)

  # Extract Fitted values: the fitted values are a list in a list. E
  fitted <- lapply(first, function(x) x[c("fitted")])
  fitted <- lapply(fitted, data.frame, stringsAsFactors = FALSE)
  fitted <- dplyr::bind_rows(fitted) # make a matrix with the fitted value of the frist stage
  # -------------------


  mydep <- paste("fitted", U, sep = "_")
  colnames(fitted) <- mydep
  second <- bind_cols(data, fitted)


  if (fex == "~0") {
    fex <- "~1"
  }
  if (method == "fe") { # no second stage FE possible if you have  method = "fe"
    b <- cbind(group, end) %>%
      group_by(group) %>%
      transmute(across(everything(), funs(. - mean(.)), .names = "{.col}dem"))
    inst_s <- paste0(paste(endog_var, collapse = "dem +"), "dem")
    second <- bind_cols(b[, -1], data, fitted)
    form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", str_sub(fen, 2, -1), "~", inst_s))
    res <- feols(form, second, cluster = group)
  } else if (method == "ols") {
    if (length(fe) == 0) {
      form <- as.Formula(paste0(".[mydep]", fex))
    } else {
      form <- as.Formula(paste0(".[mydep]", paste(fex, "|", str_sub(ffe, 2, -1))))
    }
    res <- feols(form, second, cluster = group)
  } else if (method == "2sls") {
    if (length(fe) == 0) {
      form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", str_sub(fen, 2, -1), fz))
    } else {
      form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", str_sub(ffe, 2, -1), "|", str_sub(fen, 2, -1), fz))
    }
    res <- feols(form, second, cluster = group)
  } else if (method == "be") {
    b <- cbind(group, end) %>%
      group_by(group) %>%
      transmute(across(everything(), funs(mean(.)), .names = "{.col}m"))
    second <- bind_cols(b[, -1], data, fitted)
    inst_s <- paste0(paste(endog_var, collapse = "m +"), "m")
    form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", str_sub(fen, 2, -1), "~", inst_s))
    res <- feols(form, second, cluster = group)
  } else if (method == "reoi") {
    lambda <- sapply(first, function(x) x[c("lambda_i")])
    xlist <- sapply(first, function(x) x[c("mm")])
    Xtilde <- bdiag(xlist)
    res <- list()
    for (u in U) {
      Lambda <- lapply(lambda, function(x) x[, , which(u == U)])
      # this chunk of code comes from https://github.com/cran/plm/blob/master/R/tool_ercomp.R. This uses Nerlove method and
      # with OLS using the fitted values lead to the same estiamted sigma_alpha.

      form <- as.Formula(paste0(paste0("fitted_", u), fex))
      est <- plm::plm(form, data = second, index = c("group"), model = "within")
      pdim <- plm::pdim(est)
      N <- pdim$nT$n

      s2alpha <- sum((plm::fixef(est, type = "dmean", effect = "individual"))^2 *
        pdim$Tint$Ti / pdim$nT$N) * (pdim$nT$n / (pdim$nT$n - 1))

      xlx <- Map(quad.tform, Lambda, xlist)
      omega <- Map("+", xlx, s2alpha)
      omegainv <- lapply(omega, ginv)
      Omegainv <- Matrix::bdiag(omegainv)
      Zstar <- as.matrix(Omegainv %*% model.matrix(fex, data))
      form2 <- as.Formula(paste0(form, "| Zstar-1 "))
      rr <- AER::ivreg(form2, data = second)
      rr <- coeftest(rr, vcov = vcovCL, cluster = group)
      res[[which(u == U)]] <- rr
    }
    names(res) <- mydep
  } else if (method == "gmm") {
    # https://cran.r-project.org/web/packages/momentfit/vignettes/gmmS4.pdf
    res <- list()

    for (u in U) {
      form <- as.Formula(paste0(paste0("fitted_", u), fex, "+", str_sub(fen, 2, -1)))
      model <- momentModel(form, fz, data = second, vcov = "CL", vcovOptions = list(cluster = ~group))
      rr <- summary(gmmFit(model, type = "twostep"))
      res[[which(u == U)]] <- rr
    }
    names(res) <- mydep
  } else if (method == "regmm") {
    me <- cbind(group, exo) %>%
      group_by(group) %>%
      transmute(across(everything(), funs(mean(.)), .names = "{.col}m"))

    dm <- cbind(group, exo) %>%
      group_by(group) %>%
      transmute(across(everything(), funs(. - mean(.)), .names = "{.col}dem"))

    me <- me[, -1]
    de <- dm[, -1]

    # drop zero in demeaned
    tv <- (colSums((abs(de) == 0)) == 0)
    de <- de[tv]

    Z <- cbind(de, me)

    for (u in U) {
      form <- as.Formula(paste0(paste0("fitted_", u), fex, "+", str_sub(fen, 2, -1)))
      model <- momentModel(form, Z, data = second, vcov = "CL", vcovOptions = list(cluster = ~group))
      r <- gmmFit(model, type = "twostep")
      rr <- summary(r, sandwich = TRUE, df.adj = FALSE) # whichone do we want? see 1.4.4 in https://cran.r-project.org/web/packages/momentfit/vignettes/gmmS4.pdf
      res[[which(u == U)]] <- rr
    }
    names(res) <- mydep
  } else if (method == "ht") {
    me <- cbind(group, exo) %>%
      group_by(group) %>%
      transmute(across(everything(), funs(mean(.)), .names = "{.col}m"))

    dm <- cbind(group, exo, end) %>%
      group_by(group) %>%
      transmute(across(everything(), funs(. - mean(.)), .names = "{.col}dem"))
    me <- me[, -1]
    de <- dm[, -1]

    # drop zero in demeaned
    tv <- (colSums((abs(de) == 0)) == 0)
    de <- de[tv]
    Z <- cbind(de, me)
    second <- cbind(second, Z)
    inst_s <- paste0( "~" ,paste(names(Z), collapse = "+"))
    form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", fex, "|", str_sub(fen, 2, -1), inst_s))
    res <- feols(form, second, cluster = group)
  }
  return(res)
}
