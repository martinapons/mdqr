library(purrrlyr)
library(parallel)

# rm(list = ls())

# generate some random data to try it
ddata <- NULL
x1 <- rnorm(100)
ddata <- as.data.frame(x1)
ddata$x11 <- rnorm(100)

ddata$ggroup <- rep(1:10, each = 10)
ddata$x2 <- rep(rnorm(10), each = 10)
ddata$x22 <- rep(rnorm(10), each = 10)
ddata$z <- ddata$x2 + rep(rnorm(10), each = 10)
ddata$z2 <- ddata$x22 + rep(rnorm(10), each = 10)

ddata$dep <- ddata$x1 + ddata$x2 + rnorm(100)
ddata$fe <- c(rep(1, times = 50), rep(0, times = 50))
ddata$fe2 <- c(rep(1, times = 30), rep(0, times = 70))
library(Formula)
library(formula.tools)

quantile <- c(.1, .5)

data <- ddata
mdqr <- function(formula, data, method = c("fe", "be", "reoi", "regmm", "ols", "2sls", "gmm"), quantiles = seq(0.1, 0.9, 0.1), cores = NULL) {

  formula <- Formula::as.Formula(formula)
  mm <- model.frame(as.Formula(formula(formula, lhs = 1, rhs = -5)), data)


  group <- model.frame(formula(formula, lhs = 0, rhs = 5), data)
  group_id <- names(group) # groups are define by this variable
  names(group) <- c("group")

  mm <- cbind(mm, group)

  mm %<>% as_tibble() %>% drop_na()

  # drop too small groups ----------------------------------------
  mm %<>%
    as_tibble() %>%
    group_by_at(vars(group)) %>%
    mutate(group = cur_group_id()) %>%
    ungroup()

  mm <- mm %>% add_count(group) # number of observations in each group.
  mm <- mm %>% filter(n >= 2) # remove groups with less than 25 obs

  mm %<>%
    as_tibble() %>%
    group_by_at(vars(group)) %>%
    mutate(group = cur_group_id()) %>%
    ungroup()

  G <- max(mm$group)
  # ------------------------------------------------------------

  fdep <- formula(formula, lhs = 1, rhs = 0)
  fex <- formula(formula, lhs = 0, rhs = 1)
  fen <- formula(formula, lhs = 0, rhs = 2)
  ffe <- formula(formula, lhs = 0, rhs = 4)

  y <- model.frame(fdep, mm)
  exo <- model.frame(fex, mm)
  end <- model.frame(fen, mm)
  fz <- formula(formula, lhs = 0, rhs = 3)
  z <- model.frame(fz, mm)
  fe <- model.frame(ffe, mm)

  if (dim(end)[2] > dim(z)[2] & method != "fe" & method != "ht") stop("fewer instruments than endogenous variables. If you wish to use interval instrument select method fe or ht")
  if (dim(end)[2] == 0 & dim(z)[2] > 0) stop("External instrument is specified, but there is no endogeous variable")
  if (min(tapply(y[, 1], group[, 1], var) > 0) == 0) stop("The dependent variable must vary within groups / individuals.")


  if (sum(cbind(z, group) %>% slice_rows("group") %>% dmap(var) %>% colSums() != 0) != 1) {
    stop("The instrument is not allowed to vary within groups.")
  }


  # --------------------------------------------------------------
  # Matrix of regressors (notation as in the paper X: second stage, Xtilde: first stage)
  # tv <- (cbind(exo, end, group) %>% slice_rows("group") %>% dmap(var) %>% colSums() != 0)[-1] # TRUE if time varying
  tv <-(cbind(exo, end, group) %>% slice_rows("group") %>% summarise(across(everything(), funs(n_distinct)))  %>% colSums()/G != 1)[-1]

  time_varying <- data.frame(cbind(exo, end)[, tv])
  time_varying_var <- c(names(tv)[tv == T])
  time_constant_var <- c(names(tv)[tv == F])
  dep_var <- names(y)
  fe_var <- names(fe)
  endog_var <- names(end)
  exo_var <- names(exo)
  inst_var <- names(z)
  names(time_varying) <- time_varying_var

  second <- cbind(y, exo, end, z, fe, group)


  first <- cbind(y, time_varying, group)

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
  U <- quantiles
  # first stage is the same for all estimator. Only time varying variables incldued
  # first stage estimation

  if (is.null(cores) == 1) {
    cores <- detectCores() - 1
  }
  cl <- makeCluster(cores) # Set the number of clusters
  RNGkind(kind = "L'Ecuyer-CMRG") # Set the seed for each cluster
  clusterExport(cl, c("dep_var", "time_varying_var", "md_first_stage", "datalist", "U"), environment()) # Functions needed
  clusterEvalQ(cl, {
    library(quantreg)
    library(sandwich)
    library(dtplyr)
    library(pracma)
    library(Matrix)
  })
  first <- clusterApply(cl, datalist, md_first_stage, time_varying_var, dep_var, U)
  stopCluster(cl)

  # Extract Fitted values: the fitted values are a list in a list. E
  fitted <- lapply(first, function(x) x[c("fitted")])
  fitted <- lapply(fitted, data.frame, stringsAsFactors = FALSE)
  fitted <- dplyr::bind_rows(fitted) # make a matrix with the fitted value of the frist stage
  # -------------------


  mydep <- paste("fitted", U, sep = "_")
  colnames(fitted) <- mydep
  second <- bind_cols(second, fitted)

  # pasted variables in strings:
  exo_s <- paste(exo_var, collapse = "+")
  inst_s <- paste(inst_var, collapse = "+")
  endog_s <- paste(endog_var, collapse = "+")
  fe_s <- paste(fe_var, collapse = "+")

  if (method == "fe") { # no second stage FE possible if you have  method = "fe"
    # generate internal instruments
    b <- cbind(group, end) %>%
      group_by(group) %>%
      mutate(across(everything(), funs(. - mean(.)), .names = "{.col}_dem"))

    inst_s <- paste0(paste(endog_var, collapse = "_dem +"), "_dem")
    second <- bind_cols(b, exo, fitted)

    form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~", exo_s, "|", paste0(endog_s, "~", inst_s)))
    res <- feols(form, second, cluster = group)
  } else if (method == "ols") {
    if (length(fe) == 0) {
      form <- as.Formula(paste(".[mydep]", exo_s, sep = "~"))
    } else {
      form <- as.Formula(paste0(".[mydep]~", paste(exo_s, "|", fe_s)))
    }
    res <- feols(form, second, cluster = group)
  } else if (method == "2sls") {
    if (length(fe) == 0) {
      # form <-  as.Formula(paste0(".[mydep]~", exo_s, "|", paste0(endog_s, "~", inst_s)))
      form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~", exo_s, "|", paste0(endog_s, "~", inst_s)))
    } else {
      # form <- as.Formula(paste0(".[mydep]~", exo_s, "|" , paste(fe_s, collapse = "+"), "|",paste0(endog_s, "~", inst_s) ))
      form <- as.Formula(paste0("c(", paste(mydep, collapse = ","), ")", "~", exo_s, "|", paste(fe_s, collapse = "+"), "|", paste0(endog_s, "~", inst_s)))
      # form <-  as.Formula(paste0(paste(mydep, collapse = "|"), "~", exo_s, "|", paste(fe_s, collapse = "+"), "|", paste0("(",paste(endog_var, collapse = "|"), "~", inst_s, ")"), "|", group)) for felm()
    }
    res <- feols(form, second, cluster = group)
  }
  return(res)
  # with instruments feols can't use multiple dependent variables. Thus I use felm()
}
mdqr(dep ~ x1 + x11 | x2 + x22 | z + z2 | fe + fe2 | ggroup, ddata, method = c("ols"), quantiles = c(.1, .5))


formula <- as.Formula(dep ~ x1 + x11 | x2 + x22 | z + z2 | fe + fe2 | ggroup)
