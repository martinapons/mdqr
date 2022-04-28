library(haven)
final_small <- read_dta("C:/Users/mapo/Dropbox/Uni Bern/Master/FS 2019/Workshop Econometrics II/Data/Final Data/final_small.dta")

d <- final_small
d %<>% as_tibble() %>% drop_na()

d$mygroup <- factor(with(d, interaction(month,  county_id, year)))


d %<>%
  as_tibble()%>%  group_by_at(vars(mygroup)) %>%
  mutate(group = cur_group_id()) %>%
  ungroup()


d <- d %>% add_count(group) # number of observations in each group.
d <- d %>% filter(n< 200) # remove groups with less than 25 obs
d <- d %>% filter(group < 46) # remove groups with less than 25 obs

d %<>%
  as_tibble() %>%
  group_by_at(vars(group)) %>%
  mutate(group = cur_group_id()) %>%
  ungroup()


d$mrace <- factor(d$mrace)

mdqr(bweight ~ mage+ cigs + PM25 + mrace  | 0 | 0 | county_id + month | group, data = d, method = c( "ols"), quantiles = seq(0.1, 0.9, 0.2))

  formula <- as.Formula(bweight ~ PM25| 0 | 0 | county_id + month | group)
  formula <- as.Formula(bweight ~log(mage)+ cigs + log(PM25) + factormrace | 0 | 0 | county_id + month | group)

  formula <- as.Formula(log(bweight) ~ mage * cigs + factor(mrace) + I(mage^2) + precipitation*mage |PM25 | PM10 | year:month | group)


  formula <- as.Formula(log(bweight) ~ mage +I(mage^2)+ cigs + log(PM25):precipitation +   factor(mrace) | 0 | 0 | county_id + month | group)
  mdqr(formula, data = d, method = c( "ols"), quantiles = seq(0.1, 0.9, 0.2))
  aa <- mdqr(formula, data = d, method = c( "fe"), quantiles = seq(0.2, 0.6, 0.2))
  aa <- mdqr(formula, data = d, method = c( "2sls"), quantiles = seq(0.2, 0.6, 0.2))

  mm <- model.frame(as.Formula(formula(formula, lhs = 1, rhs = -5)), data)

  formula <- as.Formula(bweight ~mage*cigs)
  lm(formula, data = d)

  # ----------------------------------
  lambda_g <- array(0, dim = c(dim(mm)[2], dim(mm)[2], length(U))) # k = number of individual level regressors
  if (is.null(c(dim(mm)[2]))){
    lambda_g <- array(0, c(1,1, length(U)))
  }



  y <- model.frame(fdep, data1)[, 1]
  for (u in U) { # For all quantiles
    qreg <- rq(y ~ mm - 1, data = data1, tau = u)
    fitted[, 1, which(U == u)] <- fitted(qreg) # Save the fitted values
    lambda_g[, , which(U == u)] <- summary.rq(qreg, se = "ker", covariance = T)$cov / dim(data1)[1]
  }
  store <- list(fitted, lambda_g, mm)

  # ----------------------------------



  # Tests

  profvis({
    formula <- as.Formula(log(bweight) ~ log(mage)  | PM25 | PM10 | county_id + month | group)
    mdqr(formula, data = d, method = c( "2sls"), quantiles = seq(0.1, 0.9, 0.2))
  })




  # FE
  formula <- as.Formula(log(bweight) ~ 0  | mage | 0 | 0 | group)
  mdqr(formula, data = d, method = c( "fe"), quantiles = seq(0.1, 0.9, 0.2))
  # OLS
  formula <- as.Formula(log(bweight) ~ mage*factor(mrace) + PM25  | 0 | 0 | county_id + month | group)
  mdqr(formula, data = d, method = c( "ols"), quantiles = seq(0.1, 0.9, 0.2))
    #BE
  formula <- as.Formula(log(bweight) ~  PM25  | mage | 0 | 0 | group)
  mdqr(formula, data = d, method = c( "be"), quantiles = seq(0.1, 0.9, 0.2))
  #RE OI
  formula <- as.Formula(log(bweight) ~  PM25 + mage | 0 | 0 | 0 | group)
  mdqr(formula, data = d, method = c( "reoi"), quantiles = seq(0.1, 0.9, 0.2))


  #GMM
  formula <- as.Formula(log(bweight) ~   0 | PM25 |CO + PM10 | 0 | group)
  mdqr(formula, data = d, method = c( "gmm"), quantiles = seq(0.1, 0.9, 0.2), clustervar = "month")

  # REGMM
  formula <- as.Formula(log(bweight) ~    PM25 + mage | 0 |0| 0 | group)
  mdqr(formula, data = d, method = c( "regmm"), quantiles = seq(0.1, 0.9, 0.2))

  # HT
  formula <- as.Formula(log(bweight) ~    PM25 + mage | PM10 + educ |0| 0 | group)
  mdqr(formula, data = d, method = c( "ht"), quantiles = seq(0.1, 0.9, 0.2))


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
  mdqr(dep ~ x1 + x11 | x2 + x22 | z + z2 | fe + fe2 | ggroup, ddata, method = c("ols"), quantiles = c(.1, .5))


  formula <- as.Formula(dep ~ x1 + x11 | x2 + x22 | z + z2 | fe + fe2 | ggroup)
