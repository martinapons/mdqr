# md_package

This package performs the Minimum Distance Estimator suggested in Melly and Pons (2023). 


`mdqr` computes the minimum distance quantile regression (MDQR) estimator suggested in Melly and Pons (2022). This estimator can be applied whenever the data has two dimensions. For instance, we may follow the same individuals during several periods (classical panel data). Or we may have individual-level data, but the treatment varies at the county level (grouped data). In both cases, mdqr estimates the effect of the regressors on the conditional quantiles of the dependent variable.
 In this help file, we focus on grouped data applications and use the corresponding terminology.

The estimator is implemented in two steps. The first stage consists of group-level quantile regressions using individual-level covariates at each quantile of interest. In the second stage, the first-stage fitted values are regressed on individual-level and group-level variables. OLS, 2SLS, GMM, FE, RE and other estimators can be used in the second stage.

In the first step, the dependent variable is regressed on all individual-level variables, i.e., on all variables that are not constant within groups.
The command automatically finds if a variable is constant or not within groups such that the user does not need to provide this information. Quantile regression is used to estimate the first stage. Therefore, it is assumed that all the variation between individuals within groups is exogenous. In principle, it would be possible to implement an estimator with instrumental variable quantile regression in the first stage, but this would require finding individual-level instrumental variables, and it would be computationally very costly.

Since the first step consists of group-level regressions, we can keep only groups where there are at least as many individuals as individual-level regressors. By default, groups with at least one remaining degree are kept. It is possible to change the cut-off with the option `n_small`.

The second step consists of a regression of the fitted from the first step on all the variables. The second stage estimator may be selected with the option method (default is `method = "ols"`). 

By construction, the sampling error of the fitted values from the first stage, which are the dependent values of the second stage, are correlated within groups. Therefore, we must cluster the standard errors at least at the level of the groups. By default, mdqr will provide standard errors clustered at this level. The option cluster allows clustering at a higher level than the groups but not at a lower level.

For the moment, the covariance between quantiles has not been implemented. It implies that the user should refrain from using these results to test hypotheses that concern several quantiles.

The most computationally intensive step is the first step. The same first step estimates can be used for several specifications in the second step. For this reason, the fitted values are saved, and if the fitted values are provided, mdqr estimates only the second stage. The fitted values can be specified using the option `fitted_values`.

The computationally intensive first-stage quantile regressions can easily be computed in parallel on several cores using the option cores. 

## Examples

### No instrument no second stage fixed effects

``` R
# Generate artificial data.
set.seed(1234)
# Generate 100 groups with 30 individuals each.
G <- 100
Ng <- 30
group <- rep(1:G, each = Ng)
# Generate a binary treatment variable.
treatment <- rep(rbinom(G, 1,0.5), each = Ng)
# Generate a group-level variable.
g_var <- rep(rnorm(G), each = Ng)
# Generate an individual-level variable.
i_var <- rnorm(Ng*G)
# Generate the outcome.
y <- 1 + treatment + g_var + i_var + rnorm(Ng * G) * (1 + 0.2 * treatment)
# Estimate the quantile treatment effect at 19 quantiles and plot the treatment effects
fit <- mdqr(y~treatment+i_var+g_var | 0 | 0 | 0 | group, method = "ols",
 data = data.frame(y, i_var, g_var, treatment, group), quantiles = seq(0.05, 0.95, 0.05), cores = 1)
fit

#>  Standard-errors: Clustered (cluster) 
#>  Dep. var.: fitted_0.05
#>              Estimate Std. Error   t value   Pr(>|t|)    
#> (Intercept) -0.445696   0.042776 -10.41922  < 2.2e-16 ***
#> treatment    0.565012   0.063361   8.91731 2.5065e-14 ***
#> i_var        1.024562   0.036588  28.00269  < 2.2e-16 ***
#> g_var        0.968654   0.031599  30.65491  < 2.2e-16 ***
#> ---
#> Dep. var.: fitted_0.1
#>              Estimate Std. Error  t value   Pr(>|t|)    
#> (Intercept) -0.229365   0.035995 -6.37215 5.9252e-09 ***
#> treatment    0.699459   0.063277 11.05393  < 2.2e-16 ***
#> i_var        1.007927   0.032191 31.31106  < 2.2e-16 ***
#> g_var        0.973306   0.031668 30.73510  < 2.2e-16 ***
#> ---
#> Dep. var.: fitted_0.15
#>              Estimate Std. Error   t value  Pr(>|t|)    
#> (Intercept) -0.005011   0.035056 -0.142932   0.88663    
#> treatment    0.694302   0.061032 11.375991 < 2.2e-16 ***
#> i_var        1.013745   0.030284 33.474860 < 2.2e-16 ***
#> g_var        0.989829   0.030706 32.236001 < 2.2e-16 ***
#> [...]

# print the summary
summary_mdqr(fit, "treatment")


#>      Quantiles  Estimate Std. Error   t value     Pr(>|t|)
#> [1,]      0.05 0.5650124 0.06336131  8.917309 2.506499e-14
#> [2,]      0.10 0.6994586 0.06327692 11.053930 5.579360e-19
#> [3,]      0.15 0.6943020 0.06103222 11.375991 1.123082e-19
#> [4,]      0.20 0.7453594 0.05077034 14.681000 1.357370e-26
#> [5,]      0.25 0.7890144 0.04914733 16.054065 2.599708e-29
#> [6,]      0.30 0.8740479 0.05289722 16.523515 3.229964e-30
#> [7,]      0.35 0.9391670 0.04751630 19.765153 3.846744e-36
#> [8,]      0.40 0.9651943 0.04876691 19.791992 3.454742e-36
#> [9,]      0.45 0.9915907 0.05025849 19.729815 4.432088e-36
#>[10,]      0.50 1.0222010 0.05010512 20.401130 3.084810e-37
#>[11,]      0.55 1.0165268 0.05273802 19.275030 2.781002e-35
#>[12,]      0.60 1.0174902 0.05411936 18.800855 1.939070e-34
#>[13,]      0.65 1.0450784 0.05787752 18.056723 4.322256e-33
#>[14,]      0.70 1.1050574 0.05376186 20.554673 1.689724e-37
#>[15,]      0.75 1.1050909 0.05644883 19.576859 8.196547e-36
#>[16,]      0.80 1.1337024 0.05388821 21.038042 2.586919e-38
#>[17,]      0.85 1.1163927 0.05824553 19.167010 4.317846e-35
#>[18,]      0.90 1.2509481 0.06387163 19.585347 7.920994e-36
#>[19,]      0.95 1.3096944 0.09438966 13.875401 5.933288e-25

# Plot the results
plot_mdqr(fit, "treatment")

```
<img width="733" alt="image" src="https://user-images.githubusercontent.com/69034843/216453849-c87bdb89-1d84-45d1-83fb-3d4f09bb4407.png">

### Grouped data, no instrument, many fixed effects in the second stage

``` R

set.seed(1234)

# Generate data for 50 states, 20 years where groups are state x years cells.
# Each group has 30 observations.
S <- 50
Y <- 20
G <- S * Y
Ng <- 30
# Generate group identifier
group <- rep(1:G, each = Ng)
# Generate 50 states (groups) with 30 years of data.
state <- rep(1:S, each = Ng*Y)
state_fe <- rep(0.2 * rnorm(S), each = Ng*Y)
first_treated <- rep(round(runif(S) * 20), each = Ng*Y)
year <- rep(rep(1:Y, each = Ng), S)
year_fe <- year * 0.001
treated <- as.numeric(year >= first_treated)
state_char <- rep(rnorm(S*Y), each = Ng)
ind_char <- rnorm(Ng * G)
# Generate the outcome. Note that we have state and year fixed effects.
y <- 1 + treated + 0.5* state_char + 0.5 * ind_char + state_fe + year_fe +
rnorm(Ng*G)*(1 + treated * 0.2)
# Generate a data frame
dat <- data.frame(y, state, state_fe, year, year_fe, treated, state_char, ind_char, group)
# Estimate the treatment effects and cluster the standard errors at the state level
fit <- mdqr(y~treated + state_char + ind_char  | 0 | 0 | state_fe + year_fe | group, method = "ols",
 data = dat, clustervar = "state", cores = 1)
fit 

#>$results
#>Standard-errors: Clustered (state) 
#>Dep. var.: fitted_0.1
#>           Estimate Std. Error t value  Pr(>|t|)    
#>treated    0.767241   0.039365 19.4906 < 2.2e-16 ***
#>state_char 0.512136   0.010100 50.7057 < 2.2e-16 ***
#>ind_char   0.511834   0.012685 40.3503 < 2.2e-16 ***
#>---
#>[...]

 
# Print Summary
summary_mdqr(fit, "treated")

#>      Quantiles  Estimate Std. Error  t value     Pr(>|t|)
#> [1,]       0.1 0.7672407 0.03936465 19.49060 9.939590e-25
#> [2,]       0.2 0.8676201 0.03512307 24.70229 2.597472e-29
#> [3,]       0.3 0.9560468 0.03047854 31.36787 4.262718e-34
#> [4,]       0.4 0.9701119 0.02992535 32.41773 9.134545e-35
#> [5,]       0.5 1.0219111 0.02948017 34.66435 3.919339e-36
#> [6,]       0.6 1.0729570 0.02943153 36.45604 3.631132e-37
#> [7,]       0.7 1.1197084 0.02723634 41.11083 1.209034e-39
#> [8,]       0.8 1.1996580 0.02844745 42.17101 3.590926e-40
#> [9,]       0.9 1.3060507 0.03166439 41.24667 1.033160e-39
 
# Plot Results
plot_mdqr(fit, "treated")


```
<img width="710" alt="image" src="https://user-images.githubusercontent.com/69034843/216457790-d98d62ad-0865-475a-b17b-baedc5d4129e.png">



## References
Melly, B., & Pons, M. (2022). Minimum Distance Estimation of Quantile Panel Data Models. 
