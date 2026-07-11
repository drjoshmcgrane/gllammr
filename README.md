# gllammr: Generalized Linear Latent and Mixed Models in R

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/gh/drjoshmcgrane/gllammr/branch/main/graph/badge.svg)](https://codecov.io/gh/drjoshmcgrane/gllammr)

An R implementation of **Generalized Linear Latent and Mixed Models (GLLAMM)** following the Stata GLLAMM framework of Rabe-Hesketh, Skrondal, and Pickles, with a fast TMB (Template Model Builder) backend.

## Overview

gllammr fits a wide class of multilevel latent variable models through one interface:

- **Multilevel GLMMs** (Gaussian, binomial, Poisson, gamma; logit/probit/cloglog links), with random slopes, nested and crossed terms, adaptive quadrature, and level-specific survey weights
- **Ordinal and multinomial models** (cumulative, adjacent-category, forward/backward continuation-ratio, partial proportional odds), with crossed/multi-term random effects
- **Item response theory**: Rasch, 2PL, 3PL (with item-specific guessing), GRM, PCM, GPCM, NRM; EM or Laplace estimation; multilevel IRT with school/cluster random effects and person-level random slopes
- **Explanatory IRT (EIRT)**: regress item difficulty, discrimination, or polytomous step thresholds on item properties (LLTM, LLTM-plus-error, latent regression, and the polytomous LPCM framework of Kim & Wilson), plus genuinely step-level covariates that vary within items across steps; single- or multilevel
- **Differential item functioning**: logistic-regression DIF with multiple grouping factors, interactions, iterative purification, and effect-size classification; and model-based IRT-likelihood-ratio DIF with anchor items, latent impact regression, and Wald or LR tests
- **Latent class analysis**: binary, categorical, and continuous indicators; order-restricted and partially ordered (poset) classes; Rasch-structured located classes (LCR); and a latent-structure comparison spanning unrestricted, monotone, invariant-item-ordering, double-monotone, LCR, and Rasch models
- **Cognitive diagnosis models**: DINA, DINO, and G-DINA with Q-matrices, attribute hierarchies, and monotonicity constraints
- **Structural equation models**: CFA, structural regressions, MIMIC, FIML for missing data, lavaan-matching fit indices (chi-square, CFI, TLI, RMSEA with CI, SRMR), standardized solutions
- **Joint mixed-response models** sharing a random effect across Gaussian/binomial/Poisson outcomes
- **Parametric survival/frailty models** (exponential, Weibull AFT)
- **Rank-ordered (exploded) logit** (Plackett-Luce) with taste-shifter random effects and partial rankings
- **NPML**: nonparametric maximum likelihood with estimated mass points replacing the normal latent distribution
- **Model comparison and inference**: `compare_models()` (AIC/BIC deltas, Akaike weights, for any mix of model classes), cluster-robust sandwich covariances, marginal (population-averaged) predictions, parametric-bootstrap `simulate()` for every class

All models are estimated by maximum likelihood via TMB automatic differentiation with the Laplace approximation (optionally adaptive Gauss-Hermite quadrature). The C++ templates compile once at install time into a single shared library — no compilation happens at run time.

## Installation

```r
# install.packages("remotes")
remotes::install_github("drjoshmcgrane/gllammr")
```

Requires R >= 4.0.0, TMB >= 1.9.0, and a C++ compiler at install time (standard for source packages; Rtools on Windows, Xcode CLT on macOS).

## Quick Start

```r
library(gllammr)

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

Note: gllammr exports its own `binomial()` (adding probit and cloglog links), which routes to a dedicated single-random-term fitter. For crossed/nested random effects, `integration = aghq()`, or level-specific weights with binary outcomes, use `stats::binomial()` to reach the general estimator.

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

A plain numeric `weights` vector is treated as observation-level frequency/probability weights. IRT, EIRT, ordinal, LCA, CDM, and survival models accept weights too.

Weighted fits reproduce duplicated-data fits exactly: EM-based fitters weight each person's log marginal likelihood, Laplace-based fitters implement integer frequency weights by exact replication, and `aghq()` weights each group's log marginal likelihood directly (supporting arbitrary non-integer level-2 weights).

### Item response theory, explanatory IRT, and DIF

```r
# Descriptive IRT (EM by default; Laplace for multilevel or SEs)
fit <- fit_irt(resp_matrix, model = "GPCM")
fit_ml <- fit_irt(resp_matrix, model = "Rasch",
                  person_data = data.frame(school = school),
                  random = ~ (1 | school))          # multilevel IRT

# Explanatory IRT: regress item parameters on item properties
fit_lltm <- fit_eirt(resp_matrix, item_data,
                     difficulty_formula = ~ complexity + word_count,
                     item_residuals = TRUE)          # LLTM + error
fit_poly <- fit_eirt(poly_matrix, item_data,
                     difficulty_formula = ~ domain,
                     threshold_formula = ~ step_type, # LPCM step regression
                     model = "PCM")
fit_step <- fit_eirt(poly_matrix, item_data,
                     difficulty_formula = ~ domain,
                     step_formula = ~ skill_demand,   # varies within item
                     step_data = step_properties,     # one row per item-step
                     model = "PCM")

# DIF: logistic-regression DIF with multiple factors and purification ...
dif <- dif_test(resp_matrix, ~ gender * language,
                person_data = person_data, purify = TRUE)
# ... or confirmatory IRT-likelihood-ratio DIF with anchor items
# (group impact handled by a latent ability regression)
dif2 <- dif_irt(resp_matrix, ~ gender, person_data = person_data,
                anchors = 1:4)
```

### Latent class, latent structure, and cognitive diagnosis

```r
# Ordered, poset, and Rasch-structured latent classes
fit_ord  <- fit_lca(Y, nclass = 3, ordering = "increasing")
fit_lcr  <- fit_lca(Y, nclass = 4, structure = "rasch")  # located classes

# Fit the full latent-structure hierarchy (UN/MON/IIO/DM/LCR/RM) at once
comp <- latent_structure_comparison(Y, nclass = 4)

# Cognitive diagnosis with a Q-matrix and attribute hierarchy
fit_cdm <- fit_cdm(Y, Q, model = "gdina",
                   hierarchy = list(c(1, 2)))   # attribute 1 before 2

# Generic model comparison for any fitted models
compare_models(fit_ord, fit_lcr, fit_cdm)
```

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

Parameter estimates, standard errors, and log-likelihoods are cross-validated against established packages — lme4, glmmTMB, ordinal, VGAM, nnet, mirt, TAM, ltm, eRm-style cross-walks, poLCA, CDM, difR, survival, npmlreg, lavaan, and clubSandwich — by an automated harness of 80+ checks covering GLMMs, ordinal/multinomial models, dichotomous and polytomous IRT, explanatory IRT (including the De Boeck & Wilson verbal-aggression benchmarks), DIF, latent class and cognitive diagnosis models, SEM (including FIML), survival, rank-ordered logit, NPML, adaptive quadrature, and sandwich standard errors:

```r
results <- gllammr_validate()   # requires the reference packages (Suggests)
```

## Documentation

```r
?gllamm
help(package = "gllammr")
browseVignettes("gllammr")
```

Thirteen vignettes: getting started, multilevel GLMMs, IRT models, multilevel IRT, explanatory IRT, DIF analysis, latent class analysis (including the latent-structure framework), cognitive diagnosis, SEM, marginal predictions, weights, advanced features, and migrating from Stata GLLAMM.

## References

- Rabe-Hesketh, S., Skrondal, A., & Pickles, A. (2004). Generalized multilevel structural equation modeling. *Psychometrika*, 69(2), 167-190.
- Skrondal, A., & Rabe-Hesketh, S. (2004). *Generalized Latent Variable Modeling: Multilevel, Longitudinal, and Structural Equation Models*. Chapman & Hall/CRC.
- De Boeck, P., & Wilson, M. (Eds.) (2004). *Explanatory Item Response Models: A Generalized Linear and Nonlinear Approach*. Springer.
- Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016). TMB: Automatic differentiation and Laplace approximation. *Journal of Statistical Software*, 70(5), 1-21.

## License

GPL-3. See [LICENSE](LICENSE).
