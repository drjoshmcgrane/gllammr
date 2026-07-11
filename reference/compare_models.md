# Compare fitted gllammr models

Builds a comparison table for any set of fitted gllammr models (or any
objects with `logLik`, `AIC`, `BIC`, and `n_obs` components):
log-likelihood, parameter count, AIC/BIC with deltas, and Akaike
weights. Works across model classes - GLMMs, IRT, latent class, CDM,
SEM, survival - because every fitter reports the same marginal
quantities.

## Usage

``` r
compare_models(..., sort_by = c("none", "AIC", "BIC"))
```

## Arguments

- ...:

  Named fitted model objects (names label the table rows; unnamed
  arguments are labelled by their expressions)

- sort_by:

  "none" (default; preserve input order), "AIC", or "BIC"

## Value

An object of class `gllammr_model_comparison`: a data frame with one row
per model.

## Details

Information criteria are only meaningful across models fitted to the
*same response data*; the function checks that `n_obs` agrees and warns
otherwise, but cannot verify that the responses themselves coincide -
that is the analyst's responsibility. For inequality-constrained models
(ordered/partially ordered classes, monotone CDMs) the parameter count
is nominal and likelihood-ratio comparisons have chi-bar-square null
distributions; treat those comparisons descriptively.

## Examples

``` r
if (FALSE) { # \dontrun{
compare_models(rasch = fit_irt(Y, "Rasch"),
               twopl = fit_irt(Y, "2PL"),
               lca3  = fit_lca(Y, nclass = 3))
} # }
```
