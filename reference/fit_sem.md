# Fit Structural Equation Models with Latent Variables

Fits a SEM with continuous indicators: a measurement model defining each
latent variable by its indicators (marker-variable identification: the
first indicator's loading is fixed at 1), an optional recursive
structural model of regressions among latent variables and on observed
covariates (MIMIC models), freely correlated exogenous latent variables,
full-information maximum likelihood for missing data, and the standard
fit-index battery.

## Usage

``` r
fit_sem(
  measurement,
  structural = NULL,
  data,
  method = c("ml", "laplace"),
  missing = c("listwise", "fiml"),
  se = TRUE,
  start = NULL,
  control = list()
)
```

## Arguments

- measurement:

  Named list of one-sided formulas defining each latent variable, e.g.
  `list(ability = ~ y1 + y2 + y3, motivation = ~ y4 + y5 + y6)`.

- structural:

  Optional list of formulas regressing latent variables on other latent
  variables and/or observed covariates, e.g.
  `list(motivation ~ ability + ses)`. The latent part must be recursive
  (no cycles); observed covariates are carried as perfectly-measured
  exogenous variables (the joint-normal formulation, equivalent to
  `lavaan` with `fixed.x = FALSE`).

- data:

  Data frame with the indicator (and covariate) variables

- method:

  Estimation method: "ml" (default; Wishart maximum likelihood on the
  sample covariance for complete data, casewise FIML under
  `missing = "fiml"`) or "laplace" (legacy full-data TMB path;
  latent-only structural models, uncorrelated exogenous factors,
  complete data).

- missing:

  "listwise" (default) or "fiml" (full-information ML over all observed
  values, assuming MAR)

- se:

  Compute standard errors (default TRUE)

- start:

  Optional starting values (laplace method only)

- control:

  Optimization control list

## Value

An object of class `gllamm_sem`. Key components: `param_table`
(estimates, SEs, z, p), `fit_measures` (chisq, df, CFI, TLI, RMSEA with
90% CI, SRMR), `latent_covariance` (exogenous block free; disturbances
diagonal), `standardized` (std.all solution), `loadings`, `structural`,
`factor_scores`, `logLik`/`AIC`/ `BIC`.

## Details

Exogenous latent variables (those with no incoming structural paths)
covary freely, as in lavaan. Likelihood-equivalent to lavaan with
`fixed.x = FALSE` when covariates are present (covariate moments are
modeled jointly), and to default lavaan otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_sem(
  measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
  structural = list(f2 ~ f1 + w),
  data = d, missing = "fiml")
summary(fit)
} # }
```
