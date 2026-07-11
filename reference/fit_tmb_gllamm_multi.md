# GLMM engine for multiple random-effects terms (crossed/nested)

Builds the lme4-style combined sparse Z mapping the full term-major
random-effects vector to observations, and fits via the glmm_multi TMB
model.

## Usage

``` r
fit_tmb_gllamm_multi(
  model_data,
  family,
  random_terms,
  start_params = NULL,
  control = list(),
  weights = NULL
)
```
