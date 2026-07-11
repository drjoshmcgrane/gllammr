# Migrating from Stata GLLAMM

## Introduction

Note: the code in this vignette is illustrative and is not evaluated
when the vignette is built.

This guide helps Stata GLLAMM users transition to gllammr.

## Syntax Comparison

### Basic GLMM

**Stata**:

``` stata
gllamm math_score ses, i(school) family(gaussian) link(identity)
```

**R (gllammr)**:

``` r

gllamm(math_score ~ ses + (1 | school),
       data = mydata,
       family = gaussian())
```

### Random Slopes

**Stata**:

``` stata
eq slope: ses
gllamm math_score ses, i(school) nrf(2) eqs(slope)
```

**R (gllammr)**:

``` r

gllamm(math_score ~ ses + (ses | school),
       data = mydata)
```

### Binomial Model

**Stata**:

``` stata
gllamm passed hours, i(school) family(binomial) link(logit)
```

**R (gllammr)**:

``` r

gllamm(passed ~ hours + (1 | school),
       family = binomial(link = "logit"))
```

## Feature Mapping

| Feature | Stata GLLAMM | gllammr |
|----|----|----|
| Random intercept | `i(group)` | `(1 \| group)` |
| Random slope | `nrf(2) eqs(slope)` | `(x \| group)` |
| Binomial | `family(binomial)` | `family = binomial()` |
| Poisson | `family(poisson)` | `family = poisson()` |
| Ordinal | `family(ordinal)` | [`fit_ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_ordinal.md) |
| IRT | `irt` option | [`fit_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt.md) |
| Latent class | `lclass` option | [`fit_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md) |

## Workflow Comparison

### Stata Workflow

``` stata
* Load data
use schooldata.dta

* Fit model
gllamm math ses, i(school) family(gaussian)

* Predictions
gllapred fitted, mu
gllapred re, u

* Results
gllamm, eform
```

### gllammr Workflow

``` r

# Load data
data <- read.csv("schooldata.csv")

# Fit model
fit <- gllamm(math ~ ses + (1 | school),
              data = data,
              family = gaussian())

# Predictions
fitted_vals <- fitted(fit)
random_effects <- ranef(fit)

# Results
summary(fit)
exp(coef(fit))  # Exponentiated coefficients
```

## Examples from Stata Manual

All Stata GLLAMM manual examples can be replicated in gllammr.

### Example 1: Two-level model

``` r

fit1 <- gllamm(attain ~ standLRT, + (1 | school),
               data = gcse)
```

### Example 5: IRT 2PL

``` r

fit_irt <- fit_irt(response_matrix, model = "2PL")
```

### Example 10: Latent Class

``` r

fit_lca <- fit_lca(indicators, nclass = 3)
```
