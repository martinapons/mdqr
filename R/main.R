library(purrrlyr)
library(parallel)
library(stringr)
library(Formula)
library(formula.tools)

quantile <- c(.1, .5)

data <- d
mdqr <- function(formula, data, method = c("fe", "be", "reoi", "regmm", "ols", "2sls", "gmm"), quantiles = seq(0.1, 0.9, 0.1), cores = NULL) {
  formula <- Formula::as.Formula(formula)

  myvar <- all.vars(formula)
  data <- select(data, all_of(myvar)) # eventually need to add cluster var.
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

  # dep_var <- str_sub(as.character(fdep), 1,  -5)

  # strsplit(as.character(ffe), " + ")

  y <- model.frame(fdep, data)
  exo <- model.frame(fex, data)
  end <- model.frame(fen, data)
  z <- model.frame(fz, data)
  fe <- model.frame(ffe, data)
  # -----------------
  dep_var <- names(y)
  fe_var <- names(fe)
  endog_var <- names(end)
  exo_var <- names(exo)
  inst_var <- names(z)


  if (dim(end)[2] > dim(z)[2] & method != "fe" & method != "ht") stop("fewer instruments than endogenous variables. If you wish to use interval instrument select method fe or ht")
  if (dim(end)[2] == 0 & dim(z)[2] > 0) stop("External instrument is specified, but there is no endogeous variable")
  if (min(tapply(y[, 1], group[, 1], var) > 0) == 0) stop("The dependent variable must vary within groups / individuals.")
  if (dim(end)[2] == 0 & method == "fe") stop("FE is used, but no endogenous variable specified.")


  if (sum(cbind(z, group) %>% slice_rows("group") %>% dmap(var) %>% colSums() != 0) != 1) {
    stop("The instrument is not allowed to vary within groups.")
  }
  # --------------------------------------------------------------
  tv <- (cbind(group, exo, end) %>% slice_rows("group") %>% summarise(across(everything(), funs(n_distinct))) %>% colSums() / G != 1)[-1]

  time_varying <- data.frame(cbind(exo, end)[, tv])
  time_varying_var <- c(names(cbind(exo, end))[tv == T]) # I don't like this two lines. I would like to take the names from tv (but they change)
  time_constant_var <- c(names(cbind(exo, end))[tv == F])

  names(time_varying) <- time_varying_var

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
  second <- bind_cols(data, fitted)


  if (fex == "~0") {
    fex <- "~1"
  }
  if (method == "fe") { # no second stage FE possible if you have  method = "fe"
    # generate internal instruments
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
  }
  return(res)
}
