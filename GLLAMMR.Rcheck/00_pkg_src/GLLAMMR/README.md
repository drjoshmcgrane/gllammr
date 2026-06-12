# GLLAMMR: Generalized Linear Latent and Mixed Models in R

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

An R implementation of **Generalized Linear Latent and Mixed Models (GLLAMM)** following the Stata GLLAMM framework of Rabe-Hesketh, Skrondal, and Pickles, with a fast TMB (Template Model Builder) backend.

## Overview

GLLAMMR fits a wide class of multilevel latent variable models through one interface:

- **Multilevel GLMMs** (Gaussian, binomial, Poisson, gamma; logit/probit/cloglog links)
- **Ordinal and multinomial models** (cumulative, adjacent-category, continuation-ratio, partial proportional odds)
- **Item response theory** (Rasch, 2PL, 3PL, GRM, PCM, GPCM, NRM), DIF analysis, explanatory IRT with item/person covariates, multilevel IRT
- **Latent class analysis**, including mixed binary/categorical/continuous indicators
- **Structural equation models** and **joint mixed-response models**
- **Parametric survival/frailty models** (exponential, Weibull)
- **Rank-ordered (exploded) logit** and **NPML mass-point** random effects

All models are estimated by maximum likelihood via TMB automatic differentiation with the Laplace approximation (optionally adaptive Gauss-Hermite quadrature). The C++ templates compile once at install time into a single shared library — no compilation happens at run time.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yourusername/GLLAMMR")
```

Requires R >= 4.0.0, TMB >= 1.9.0, and a C++ compiler at install time (standard for source packages; Rtools on Windows, Xcode CLT on macOS).

## Quick Start

```r
library(GLLAMMR)

# Gaussian GLMM with random intercept
fit1 <- gllamm(y ~ x + (1 | group), data = mydata, family = gaussian())

# Logistic GLMM
fit2 <- gllamm(y ~ x + (1 | group), data = mydata, family = binomial())

# Random intercept and slope
fit3 <- gllamm(y ~ x + (x | group), data = mydata)

summary(fit1)
fixef(fit1)     # fixed effects
ranef(fit1)     # random effects (empirical Bayes)
VarCorr(fit1)   # variance components
```

### Random slopes, nested, and crossed random effects

```r
# Random intercept + slope with full covariance (sleepstudy example)
data(sleepstudy, package = "lme4")
fit <- gllamm(Reaction ~ Days + (Days | Subject), data = sleepstudy)

# Nested: classes within schools
fit_nested <- gllamm(score ~ ses + (1 | school/class), data = school_data)

# Crossed: persons crossed with items
fit_crossed <- gllamm(correct ~ x + (1 | person) + (1 | item),
                      data = test_data, family = stats::binomial())
```

| Syntax | Description |
|--------|-------------|
| `(1 \| g)` | Random intercept |
| `(x \| g)` | Correlated random intercept and slope |
| `(x \|\| g)` | Uncorrelated random intercept and slope |
| `(1 \| g1/g2)` | Nested random effects |
| `(1 \| g1) + (1 \| g2)` | Crossed random effects |

Note: GLLAMMR exports its own `binomial()` (adding probit and cloglog links), which routes to a dedicated single-random-term fitter. For crossed/nested random effects, `integration = aghq()`, or level-specific weights with binary outcomes, use `stats::binomial()` to reach the general estimator.

### Unified family dispatch

Latent variable models are available through `gllamm()` with special family objects, or directly through `fit_*()` functions:

```r
# IRT: first argument is the persons x items response matrix
fit_2pl <- gllamm(resp_matrix, family = irt("2PL"))
fit_2pl <- fit_irt(resp_matrix, model = "2PL")        # equivalent
fit_2pl$item_parameters
abilities(fit_2pl)

# Latent class analysis
fit_lca <- gllamm(indicator_matrix, family = lca(nclass = 3))
fit_lca$class_probs
fit_lca$posterior

# Ordinal and multinomial
fit_ord <- gllamm(rating ~ temp + (1 | judge), data = wine,
                  family = ordinal("logit"))
fit_mnl <- gllamm(choice ~ x + (1 | region), data = d, family = multinomial())
```

### Adaptive quadrature

The default Laplace approximation can be replaced with adaptive Gauss-Hermite quadrature for small-cluster binary/count data:

```r
fit <- gllamm(y ~ x + (1 | g), data = d, family = stats::binomial(),
              integration = aghq(15))
```

### Survey weights

Level-specific (pseudo-likelihood) weights, as in Stata GLLAMM's `pweight()`:

```r
fit <- gllamm(y ~ x + (1 | cluster), data = d, family = stats::binomial(),
              weights = list(level1 = obs_weights, level2 = cluster_weights))
```

A plain numeric `weights` vector is treated as observation-level frequency/probability weights. IRT, ordinal, LCA, and survival models accept weights too.

### Other model classes

```r
# Parametric survival with shared frailty
fit_surv <- fit_survival(Surv(time, status) ~ x + (1 | center),
                         data = d, distribution = "weibull")

# Joint mixed responses sharing a random effect
fit_mr <- fit_mixed(formulas = list(gaussian = yc ~ x, binomial = yb ~ x),
                    random = ~ (1 | grp), data = d)

# Structural equation model
fit_sem <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                      f2 = ~ y1 + y2 + y3),
                   structural = list(f2 ~ f1), data = d)

# Rank-ordered (exploded) logit with random coefficients
fit_rank <- fit_rank(rank ~ price + quality, case = ~ id,
                     random = ~ (0 + price | region), data = d)

# Nonparametric maximum likelihood (discrete mass points)
fit_np <- fit_npml(y ~ x + (1 | grp), data = d, k = 3, family = binomial())
```

### Robust standard errors

```r
vcov(fit, type = "sandwich")   # cluster-robust (sandwich) covariance
```

### Predictions and simulation

```r
pred_cond <- predict(fit, type = "response")             # conditional on u = 0
pred_marg <- predict(fit, type = "marginal", n_sim = 1000) # population-averaged
sims <- simulate(fit, nsim = 100)                        # data frame of draws
```

For nonlinear links, conditional and marginal predictions differ; marginal predictions integrate over the random effects by Monte Carlo (`se.fit = TRUE` available). See `vignette("marginal-predictions")`.

## Validation

Parameter estimates, standard errors, and log-likelihoods are cross-validated against established packages — lme4, glmmTMB, ordinal, mirt, poLCA, npmlreg, and lavaan — by an automated suite of 49 checks covering Gaussian/binomial/Poisson/gamma GLMMs, ordinal models, Rasch/2PL/GRM, latent class models, survival, SEM, NPML, and adaptive quadrature:

```r
results <- gllammr_validate()   # requires the reference packages (Suggests)
```

## Documentation

```r
?gllamm
help(package = "GLLAMMR")
browseVignettes("GLLAMMR")
```

Vignettes: getting started, multilevel GLMMs, IRT models, multilevel IRT, latent class analysis, marginal predictions, weights, advanced features, and migrating from Stata GLLAMM.

## References

- Rabe-Hesketh, S., Skrondal, A., & Pickles, A. (2004). Generalized multilevel structural equation modeling. *Psychometrika*, 69(2), 167-190.
- Skrondal, A., & Rabe-Hesketh, S. (2004). *Generalized Latent Variable Modeling: Multilevel, Longitudinal, and Structural Equation Models*. Chapman & Hall/CRC.
- Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016). TMB: Automatic differentiation and Laplace approximation. *Journal of Statistical Software*, 70(5), 1-21.

## License

GPL-3. See [LICENSE](LICENSE).
