# Generic Fit Statistics Function

Compute model-specific fit statistics for GLLAMM objects

## Usage

``` r
fit(object, ...)
```

## Arguments

- object:

  A fitted model object

- ...:

  Additional arguments passed to methods

## Value

An object of class `fit_statistics` with model-specific components

## Details

The `fit()` function provides comprehensive, model-specific fit
statistics:

**GLMM models:**

- Log-likelihood, AIC, BIC

- R-squared (marginal and conditional for Gaussian models)

- Intraclass correlation coefficient (ICC)

**IRT models:**

- Log-likelihood, AIC, BIC

- Item fit statistics (S-X^2)

- Person fit statistics (outfit/infit)

- Reliability estimates

- Test information function

**Latent Class Analysis:**

- Log-likelihood, AIC, BIC

- Entropy (classification quality)

- Class proportions

- Average posterior probabilities (APPA)

**Ordinal models:**

- Log-likelihood, AIC, BIC

- Pseudo-R^2 (McFadden)

- Proportional odds test (for PO/probit models)

## Examples

``` r
if (FALSE) { # \dontrun{
# GLMM
fit1 <- gllamm(y ~ x + (1 | group), data = data)
fit(fit1)

# IRT
fit2 <- fit_irt(responses, model = "2PL")
fit(fit2, compute_item_fit = TRUE)

# LCA
fit3 <- fit_lca(indicators, nclass = 3)
fit(fit3)

# Ordinal
fit4 <- fit_ordinal(rating ~ x + (1 | id), data = data)
fit(fit4, test_po = TRUE)
} # }
```
