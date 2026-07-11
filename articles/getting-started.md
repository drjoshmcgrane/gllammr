# Getting Started with gllammr

## Introduction

gllammr (Generalized Linear Latent and Mixed Models in R) provides a
comprehensive framework for fitting multilevel latent variable models.
This vignette introduces the basic functionality.

## Installation

``` r

# From CRAN (once published):
# install.packages("gllammr")
# Development version:
# remotes::install_github("drjoshmcgrane/gllammr")
library(gllammr)
#> 
#> Attaching package: 'gllammr'
#> The following object is masked from 'package:stats':
#> 
#>     binomial
```

## Basic Random Intercept Model

The simplest GLLAMM model is a random intercept model:

``` math
y_{ij} = \beta_0 + \beta_1 x_{ij} + u_j + \epsilon_{ij}
```

where: - $`y_{ij}`$ is the response for observation $`i`$ in group
$`j`$ - $`x_{ij}`$ is a predictor - $`u_j \sim N(0, \sigma_u^2)`$ is the
random intercept for group $`j`$ - $`\epsilon_{ij} \sim N(0, \sigma^2)`$
is the residual error

### Simulate Data

``` r

set.seed(123)

# Parameters
n_groups <- 30
n_per_group <- 10
beta_0 <- 2.0
beta_1 <- -0.5
sigma_u <- 0.8
sigma <- 1.2

# Generate data
group <- rep(1:n_groups, each = n_per_group)
x <- rnorm(n_groups * n_per_group)
u <- rnorm(n_groups, sd = sigma_u)
y <- beta_0 + beta_1 * x + u[group] + rnorm(n_groups * n_per_group, sd = sigma)

data <- data.frame(y = y, x = x, group = group)

# View first few rows
head(data)
#>            y           x group
#> 1  2.9291144 -0.56047565     1
#> 2  0.1167742 -0.23017749     1
#> 3 -0.2174732  1.55870831     1
#> 4  3.2156133  0.07050839     1
#> 5  1.8160280  0.12928774     1
#> 6 -1.8923936  1.71506499     1
```

### Fit Model

``` r

fit <- gllamm(y ~ x + (1 | group),
              data = data,
              family = gaussian())

# Print results
print(fit)
#> Generalized Linear Latent and Mixed Model
#> 
#> Call:
#> gllamm(formula = y ~ x + (1 | group), data = data, family = gaussian())
#> 
#> Family: gaussian 
#> Link: identity 
#> 
#> Random effects:
#>   Groups: group 
#>     Number of groups: 30 
#> 
#> Fixed effects:
#> (Intercept)           x 
#>      2.0317     -0.3906 
#> 
#> Log-likelihood: -509.68 
#> AIC: 1027.35   BIC: 1042.17
```

### Model Summary

``` r

summary(fit)
#> Generalized Linear Latent and Mixed Model
#> 
#> Call:
#> gllamm(formula = y ~ x + (1 | group), data = data, family = gaussian())
#> 
#> Family: gaussian 
#> Link: identity 
#> 
#> Fixed effects:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.03169    0.16485  12.324  < 2e-16 ***
#> x           -0.39057    0.07685  -5.082 3.73e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random effects:
#>   Groups: group 
#>     Number of groups: 30 
#>     Variance: 0.6675 
#>     Std.Dev.: 0.817 
#> 
#> Number of observations: 300 
#> Log-likelihood: -509.68 
#> AIC: 1027.35 
#> BIC: 1042.17
```

## Extracting Components

gllammr provides several functions to extract model components:

``` r

# Fixed effects coefficients
fixef(fit)
#> (Intercept)           x 
#>   2.0316884  -0.3905713

# Random effects (BLUPs)
ranef(fit)
#> [[1]]
#>          u 
#> -0.6829593 
#> 
#> [[2]]
#>          u 
#> -0.5497182 
#> 
#> [[3]]
#>          u 
#> -0.7368353 
#> 
#> [[4]]
#>          u 
#> -0.5300653 
#> 
#> [[5]]
#>          u 
#> -0.4697228 
#> 
#> [[6]]
#>         u 
#> 0.2134484 
#> 
#> [[7]]
#>         u 
#> -1.687158 
#> 
#> [[8]]
#>          u 
#> -0.1287837 
#> 
#> [[9]]
#>         u 
#> 0.6138465 
#> 
#> [[10]]
#>        u 
#> 1.751324 
#> 
#> [[11]]
#>         u 
#> 0.7865736 
#> 
#> [[12]]
#>         u 
#> 0.1870099 
#> 
#> [[13]]
#>          u 
#> -0.9071237 
#> 
#> [[14]]
#>          u 
#> 0.01169199 
#> 
#> [[15]]
#>           u 
#> -0.05680421 
#> 
#> [[16]]
#>         u 
#> 0.9013834 
#> 
#> [[17]]
#>           u 
#> -0.06871744 
#> 
#> [[18]]
#>         u 
#> -1.278806 
#> 
#> [[19]]
#>         u 
#> 0.9606721 
#> 
#> [[20]]
#>         u 
#> 0.8680124 
#> 
#> [[21]]
#>         u 
#> 0.3244595 
#> 
#> [[22]]
#>         u 
#> 0.9671547 
#> 
#> [[23]]
#>          u 
#> -0.8823734 
#> 
#> [[24]]
#>         u 
#> 0.4093617 
#> 
#> [[25]]
#>          u 
#> -0.4363732 
#> 
#> [[26]]
#>          u 
#> 0.07437115 
#> 
#> [[27]]
#>          u 
#> -0.2280285 
#> 
#> [[28]]
#>          u 
#> 0.06036149 
#> 
#> [[29]]
#>         u 
#> 0.6142597 
#> 
#> [[30]]
#>          u 
#> -0.1004529

# Variance components
VarCorr(fit)
#> Random effects variance components:
#> 
#>  Group: group 
#>        [,1]
#> [1,] 0.6675

# Fitted values
head(fitted(fit))
#>         1         2         3         4         5         6 
#> 1.5676348 1.4386298 0.7399424 1.3211906 1.2982330 0.6788739

# Residuals
head(residuals(fit))
#>          1          2          3          4          5          6 
#>  1.3614796 -1.3218557 -0.9574156  1.8944227  0.5177949 -2.5712676

# Log-likelihood
logLik(fit)
#> 'log Lik.' -509.6767 (df=4)

# AIC and BIC
AIC(fit)
#> [1] 1027.353
BIC(fit)
#> [1] 1042.169
```

## Random Intercept and Slope

For a more complex model with both random intercepts and slopes:

``` math
y_{ij} = \beta_0 + \beta_1 x_{ij} + u_{0j} + u_{1j}x_{ij} + \epsilon_{ij}
```

``` r

# Simulate data with random slopes
set.seed(456)
n_groups <- 20
n_per_group <- 15

group <- rep(1:n_groups, each = n_per_group)
x <- rnorm(n_groups * n_per_group)

# Random intercepts and slopes
u0 <- rnorm(n_groups, sd = 1.0)
u1 <- rnorm(n_groups, sd = 0.5)

y <- 3.0 + 0.8 * x + u0[group] + u1[group] * x + rnorm(n_groups * n_per_group, sd = 1.5)

data_rs <- data.frame(y = y, x = x, group = group)

# Fit model with random slope
fit_rs <- gllamm(y ~ x + (x | group), data = data_rs)

summary(fit_rs)
#> Generalized Linear Latent and Mixed Model
#> 
#> Call:
#> gllamm(formula = y ~ x + (x | group), data = data_rs)
#> 
#> Family: gaussian 
#> Link: identity 
#> 
#> Fixed effects:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   3.3644     0.1782  18.881  < 2e-16 ***
#> x             0.8784     0.1207   7.279 3.36e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random effects:
#>   Groups: group 
#>     Number of groups: 20 
#>     Variance: 0.4857 0.0477 0.0477 0.1198 
#>     Std.Dev.: 0.6969 0.2184 0.2184 0.3461 
#> 
#> Number of observations: 300 
#> Log-likelihood: -563.38 
#> AIC: 1138.75 
#> BIC: 1160.97
```

## Three-Level Model

gllammr supports models with multiple levels of nesting:

``` r

# Students nested in classes nested in schools
# y_ijk = β0 + β1*x_ijk + u_school[k] + u_class[jk] + ε_ijk

# Simulate three-level data
set.seed(789)
n_schools <- 15
n_classes_per_school <- 3
n_students_per_class <- 8

n_obs <- n_schools * n_classes_per_school * n_students_per_class

school <- rep(1:n_schools, each = n_classes_per_school * n_students_per_class)
class <- rep(1:(n_schools * n_classes_per_school), each = n_students_per_class)
x <- rnorm(n_obs)

u_school <- rnorm(n_schools, sd = 1.2)
u_class <- rnorm(n_schools * n_classes_per_school, sd = 0.8)

y <- 5.0 + 0.5 * x + u_school[school] + u_class[class] + rnorm(n_obs, sd = 1.0)

data_3level <- data.frame(y = y, x = x, school = school, class = class)

# Fit three-level model
fit_3level <- gllamm(y ~ x + (1 | school/class), data = data_3level)

summary(fit_3level)
#> Generalized Linear Latent and Mixed Model
#> 
#> Call:
#> gllamm(formula = y ~ x + (1 | school/class), data = data_3level)
#> 
#> Family: gaussian 
#> Link: identity 
#> 
#> Fixed effects:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  5.04577    0.24725  20.408   <2e-16 ***
#> x            0.47704    0.05129   9.301   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random effects:
#>   Groups: school 
#>     Number of groups: 15 
#>     Variance: 0.7007 
#>     Std.Dev.: 0.8371 
#>   Groups: school:class 
#>     Number of groups: 45 
#>     Variance: 0.5377 
#>     Std.Dev.: 0.7333 
#> 
#> Number of observations: 360 
#> Log-likelihood: -540.23 
#> AIC: 1090.45 
#> BIC: 1109.88
```

## Prediction

Generate predictions from fitted models:

``` r

# Predictions with random effects (default)
pred_full <- predict(fit)

# Population-level predictions (no random effects)
pred_pop <- predict(fit, re.form = NA)

# Compare
head(data.frame(
  observed = data$y,
  fitted_with_re = pred_full,
  fitted_no_re = pred_pop
))
#>     observed fitted_with_re fitted_no_re
#> 1  2.9291144      1.5676348     2.250594
#> 2  0.1167742      1.4386298     2.121589
#> 3 -0.2174732      0.7399424     1.422902
#> 4  3.2156133      1.3211906     2.004150
#> 5  1.8160280      1.2982330     1.981192
#> 6 -1.8923936      0.6788739     1.361833
```

## Simulation

Simulate new data from a fitted model:

``` r

# Single simulation
sim1 <- simulate(fit, nsim = 1, seed = 999)

# Multiple simulations
sim10 <- simulate(fit, nsim = 10, seed = 999)
head(sim10)
#>       sim_1      sim_2     sim_3       sim_4     sim_5    sim_6      sim_7
#> 1  4.914727  3.1421377 0.3489750 2.340022283 2.4983499 2.564305  2.8767315
#> 2  2.621796  1.9380250 2.8661383 2.476861064 1.7418699 2.107857  3.3412807
#> 3  1.410591  1.9005236 0.1088693 0.005004313 0.2123433 2.327248  0.7005704
#> 4  3.086528 -0.9005855 1.0240178 0.068504049 2.4805534 3.077065  3.9581190
#> 5  1.451191  2.7703469 0.7553424 1.139294259 1.7937410 2.725704  2.5233019
#> 6 -1.435999 -0.0404280 1.2646073 3.371964218 3.8798912 1.818402 -0.2447819
#>         sim_8       sim_9     sim_10
#> 1 -0.02631213  0.04180876  1.5782663
#> 2  0.94410549  0.94458925  0.1881216
#> 3  2.13780982  0.43387288 -0.5763878
#> 4  2.19968613  2.26827264  0.1071761
#> 5  1.93557923  0.37291674 -1.2237347
#> 6  0.44779099 -0.10227832 -1.9629194
```

## The full model space

Every model class in the package runs through
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md);
the family object selects the model:

| Model class | Family / call shape | Vignette |
|----|----|----|
| GLMMs (crossed/nested REs, slopes) | [`gaussian()`](https://rdrr.io/r/stats/family.html), [`binomial()`](https://drjoshmcgrane.github.io/gllammr/reference/binomial.md), [`poisson()`](https://rdrr.io/r/stats/family.html), [`Gamma()`](https://rdrr.io/r/stats/family.html) | `multilevel-glmm` |
| Ordinal (5 links), multinomial | `ordinal(link)`, [`multinomial()`](https://drjoshmcgrane.github.io/gllammr/reference/multinomial.md) | `multilevel-glmm` |
| IRT (Rasch … NRM) | `gllamm(resp, family = irt(model))` | `irt-models` |
| Explanatory IRT | `family = eirt(item_data, ~ ...)` | `explanatory-irt` |
| Multilevel IRT | [`irt()`](https://drjoshmcgrane.github.io/gllammr/reference/irt.md) + `random = ~ (1 | group)` | `multilevel-irt` |
| DIF (screening + confirmatory) | [`dif_test()`](https://drjoshmcgrane.github.io/gllammr/reference/dif_test.md), [`dif_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/dif_irt.md) | `dif-analysis` |
| Latent classes (+ ordered, poset) | `family = lca(nclass, ordering)` | `latent-class` |
| Cognitive diagnosis | `family = cdm(Q, model)` | `cognitive-diagnosis` |
| SEM (CFA, MIMIC, FIML) | `gllamm(df, family = sem(...))` | `sem-models` |
| Mixed responses, survival, ranks, NPML, AGHQ | [`mixed_response()`](https://drjoshmcgrane.github.io/gllammr/reference/mixed_response.md), [`survival_family()`](https://drjoshmcgrane.github.io/gllammr/reference/survival_family.md), [`ranking()`](https://drjoshmcgrane.github.io/gllammr/reference/ranking.md), `integration = npml(k)/aghq(k)` | `advanced-features` |
| Survey weights | `weights = list(level1 = , level2 = )` | `weights` |
| Marginal predictions | `predict(fit, type = "marginal")` | `marginal-predictions` |
| Stata gllamm migration | syntax crosswalk | `stata-migration` |

Cross-package validation results for all of these are generated by
[`gllammr_validate()`](https://drjoshmcgrane.github.io/gllammr/reference/gllammr_validate.md).

## References

- Rabe-Hesketh, S., Skrondal, A., & Pickles, A. (2004). GLLAMM Manual.
- Skrondal, A., & Rabe-Hesketh, S. (2004). Generalized Latent Variable
  Modeling.
