# Fit Joint Models for Mixed Response Types

Fits a joint model for up to three outcomes of different types
(gaussian, binomial, poisson) measured on the same observations and
sharing a common random effect. The shared random effect induces
association between the outcomes.

## Usage

``` r
fit_mixed(formulas, random, data, start = NULL, control = list())
```

## Arguments

- formulas:

  Named list of up to three fixed-effects formulas, with names among
  "gaussian", "binomial", "poisson", e.g.
  `list(gaussian = y1 ~ x1, binomial = y2 ~ x2)`.

- random:

  One-sided random-intercept formula, e.g. `~ (1 | group)`

- data:

  Data frame containing all outcome and covariate variables; rows must
  be complete for every supplied outcome

- start:

  Optional starting values

- control:

  Optimization control list

## Value

An object of class `gllamm_mixed`

## Details

The shared random effect enters every outcome's linear predictor with
loading 1 (a common-intercept joint model). Outcome-specific loadings
are not yet supported.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_mixed(
  formulas = list(gaussian = biomarker ~ age,
                  binomial = event ~ age + treatment),
  random = ~ (1 | patient),
  data = d)
} # }
```
