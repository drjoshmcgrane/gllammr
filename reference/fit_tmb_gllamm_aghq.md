# Fit a two-level GLMM by adaptive Gauss-Hermite quadrature

R driver for the glmm_aghq TMB objective: alternates parameter
optimization with updates of the per-group adaptation centers and scales
(posterior modes and curvatures), the classic adaptive quadrature
scheme.

## Usage

``` r
fit_tmb_gllamm_aghq(
  model_data,
  family,
  random_terms,
  k = 15,
  start_params = NULL,
  control = list(),
  weights = NULL,
  max_adapt = 5,
  adapt_tol = 1e-04
)
```
