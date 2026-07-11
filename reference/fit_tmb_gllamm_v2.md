# Enhanced interface to TMB for gllammr models

Supports random intercepts, random slopes, and multiple GLM families

## Usage

``` r
fit_tmb_gllamm_v2(
  model_data,
  family,
  random_terms,
  start_params = NULL,
  control = list(),
  weights = NULL
)
```

## Arguments

- model_data:

  List from make_model_matrices()

- family:

  GLM family object

- random_terms:

  Parsed random effects terms

- start_params:

  Optional starting values

- control:

  Control parameters for optimization

## Value

TMB fit object with parameter estimates
